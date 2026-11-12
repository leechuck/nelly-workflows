"""Bounded tool-calling LLM layer around the deterministic ACMG engine.

Design §4.3 step 4 + §8 decision 4 + §9 decision 8 implemented *in code*
(fail-closed, not by prompt):

* the model reviews tier 0-2 candidates plus the top 5 overall, at most 8
  variants, 12 tool calls and 6 min per variant, 50 min global budget;
* ``validate_submission`` enforces: codes ⊆ the 28 ACMA codes, no
  ``not_evaluable`` criterion may become ``met`` (nor a retired
  ``not_applicable`` one be resurrected), any change must carry a non-empty
  reason naming an Evidence field, strengths ∈ allowed set, rationale must
  not cite a PMID or any number absent from the evidence pack;
* rejected / failed / timed-out assessments fall back to the deterministic
  result; the narrative is written only from the FINAL criteria list and the
  classification is decided by classify.py, never by the model.

Entry points (WP4 scope note):
    run_agent(ev, det, backend, ...)       -> AgentResult for one variant
    run_all(candidates, backend, ...)      -> (candidates, global_status)
    run_llm_layer(candidates, backend, ...) -> candidates  (same dicts)
"""
from __future__ import annotations

import re
import time
from dataclasses import dataclass, field, fields
from typing import Any

from acmg_agent.llm_backends import (  # noqa: F401 (re-export Message for callers)
    Message,
    _build_trace_entry,
    append_llm_trace,
)
from acmg_agent.model import Criterion, Evidence
from acmg_agent.policy import RECOMMENDATION_RULE, find_recommendations, strip_recommendations, violation_message
from acmg_agent.tools import (
    REPEAT_READ_TOOL,
    REPEAT_SUBMIT_TOOL,
    SUBMIT_TOOL,
    EvidencePack,
    RepeatPack,
    all_tool_specs,
    build_prompt_context,
    repeat_tool_specs,
)

ALL_CRITERION_CODES = frozenset({
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
})
ALLOWED_STATUSES = frozenset({"met", "not_met", "not_evaluable", "not_applicable"})
ALLOWED_STRENGTHS = frozenset({"stand_alone", "very_strong", "strong", "moderate", "supporting", None})
EVIDENCE_FIELDS = tuple(f.name for f in fields(Evidence))

# Run-level (llm.status) vocabulary used by report.py: ok|rejected|failed|timeout|partial
OFF_STATUS = "off"
GLOBAL_STATUS_OK = "ok"
GLOBAL_STATUS_PARTIAL = "partial"
GLOBAL_STATUS_FAILED = "failed"

MAX_TOOL_CALLS = 12
PER_VARIANT_SECONDS = 360       # 6 min
GLOBAL_SECONDS = 3000           # 50 min
MAX_VARIANTS = 8

PMC_RE = re.compile(r"PMID\s?\d+|pubmed\s?\d+", re.IGNORECASE)
NUMBER_RE = re.compile(r"\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")


AGENT_SYSTEM_PROMPT = (
    "You are a clinical genetics reviewer revising a deterministic ACMG/AMP 2015 "
    "assessment produced by a rule engine. You may only adjust the engine's "
    "criteria using information from the evidence pack.\n"
    "Contract:\n"
    "1. Use the read tools to inspect the pack, then submit your assessment with "
    "the submit_assessment tool. Unchanged criteria may be omitted.\n"
    "2. You may change the strength of a criterion the engine marked met, or mark "
    "a met criterion not_met, but ONLY with a reason naming the evidence-pack "
    "field you relied on.\n"
    "3. You may never mark a not_evaluable criterion as met, never invent "
    "criterion codes, and never introduce evidence not present in the pack.\n"
    "4. Every claim you make must cite a value present in the pack. Do not cite "
    "PMIDs, literature, or any number that is not in the pack.\n"
    "5. ClinVar interpretation labels are evidence for PS1/PM5 same-codon or "
    "same-variant lookups only; they are never a classification shortcut.\n"
    "6. The final classification (pathogenic vs benign vs inconclusive) is "
    "computed by the calling code, not by you. In your rationale, describe the "
    "criteria only.\n"
    "7. Follow ClinGen SVI rules encoded by the engine: PM2/BS1/BA1 rarity uses "
    "non-bottlenecked max_af (never overturn PM2 using founder/bottlenecked "
    "max_af_any_pop such as ASJ/FIN/AMI/MID), and PM3 at supporting applies to "
    "homozygous (zygosity='hom') variants in AR genes even when other_var_in_gene "
    "is empty. Keep internal reasoning concise so tool calls finish promptly.\n"
    "If you cannot improve the engine's assessment with pack evidence, submit "
    "no changes and say so in the rationale.\n"
    "8. " + RECOMMENDATION_RULE
)


@dataclass
class AgentResult:
    """Outcome of one per-variant review loop."""
    criteria: list[Criterion]
    rationale: str = ""
    status: str = "ok"            # ok | ok_truncated | rejected | failed | timeout
    narrative: str | None = None
    adjusted: list[str] = field(default_factory=list)
    n_tool_calls: int = 0
    redacted: bool = False        # sentences removed under the recommendation policy


# --------------------------------------------------------------------------
# tool-calling helpers
# --------------------------------------------------------------------------
NARRATIVE_SYSTEM_PROMPT = (
    "Write the clinical summary paragraph for this variant assessment. Grounded "
    "rules: mention only gene names and only evidence values present in the "
    "supplied material; no literature, no PMIDs, no numbers absent from the "
    "material; 4-8 sentences; state the criteria that drive the result and name "
    "the strongest near-tie as a differential when one exists. The final "
    "classification is decided by the calling code, so phrase findings evidence "
    "in neutral terms. The criteria list below is already finalized: keep "
    "internal reasoning brief and output a complete paragraph ending with a "
    "period (.). " + RECOMMENDATION_RULE
)

TERMINAL_PUNCT_RE = re.compile(r'[.!?]["\'”’)\]]*\s*$')


def _ends_with_terminal_punctuation(text: str | None) -> bool:
    """Return True iff *text* is non-empty and ends with sentence-terminal punctuation."""
    if not text or not text.strip():
        return False
    return bool(TERMINAL_PUNCT_RE.search(text.strip()))


def _call_backend_chat(
    backend: Any,
    messages: list[Message],
    tools: list[dict],
    *,
    max_tokens: int | None = None,
    trace_meta: dict[str, Any] | None = None,
    attempt: str = "initial",
) -> Message:
    """Invoke ``backend.chat`` with optional ``max_tokens``/``trace_meta`` and ensure trace logging."""
    t0 = time.monotonic()
    try:
        reply = backend.chat(messages, tools, max_tokens=max_tokens, trace_meta=trace_meta)
    except TypeError:
        reply = backend.chat(messages, tools)
    elapsed = time.monotonic() - t0
    if not getattr(reply, "_trace_logged", False):
        eff_max_tokens = max_tokens if (max_tokens is not None and max_tokens > 0) else getattr(backend, "max_tokens", None)
        entry = _build_trace_entry(
            backend_name=str(getattr(backend, "backend_name", type(backend).__name__)),
            model_name=getattr(backend, "model_name", None),
            session_affinity=getattr(backend, "session_affinity", None),
            url=getattr(backend, "url", None),
            max_tokens=eff_max_tokens,
            n_messages=len(messages),
            n_tools=len(tools or []),
            resp=None,
            msg=reply,
            elapsed_seconds=elapsed,
            trace_meta=trace_meta,
            attempt=attempt,
        )
        append_llm_trace(entry, trace_path=getattr(backend, "trace_path", None))
    return reply


def _backfill_narrative_tool_messages(
    messages: list[Any], tool_specs: list[dict]
) -> list[Any]:
    """Reusable hook point for tests that want to trim context; pass-through."""
    return messages


# --------------------------------------------------------------------------
# contract enforcement
# --------------------------------------------------------------------------
def _extract_numbers(text: str, strip_codes: bool = False) -> set[float]:
    if strip_codes:
        for code in ALL_CRITERION_CODES:
            text = text.replace(code, " ")
    out: set[float] = set()
    for tok in NUMBER_RE.findall(text):
        try:
            out.add(float(tok))
        except ValueError:  # NUMBER_RE guarantees parseability; belt and braces
            continue
    return out


def _numbers_allowed(source: str, allowed: set[float] | None) -> bool:
    if allowed is None:
        return True
    return all(value in allowed for value in _extract_numbers(source, strip_codes=True))


def _mentions_evidence_field(reason: str) -> bool:
    low = reason.lower()
    for field_name in EVIDENCE_FIELDS:
        if re.search(rf"\b{re.escape(field_name.lower())}\b", low):
            return True
    return False


def validate_submission(
    det: list[Criterion],
    submitted: dict,
    *,
    allowed_numbers: set[float] | None = None,
) -> tuple[bool, str]:
    """Enforce the agent contract. Returns (ok, reason).

    Rules (design §4.3 step 4, §8.4, §9.8):
      * submitted must contain a criteria list and a non-empty rationale;
      * every criterion code must be one of the 28 canonical codes;
      * statuses ∈ {met, not_met, not_evaluable, not_applicable};
      * ``not_evaluable`` may never become ``met`` (and retired
        ``not_applicable`` codes may never be activated);
      * ``met`` criteria carry a strength from the allowed set, all other
        statuses carry no strength;
      * any change versus the deterministic result requires a non-empty
        reason that names a field of Evidence;
      * reasons and the rationale must not cite a PMID or any number that is
        not present in the evidence pack (``allowed_numbers=None`` disables
        the numeric check so plain validator tests can run).
    """
    if not isinstance(submitted, dict):
        return (False, "submission must be an object")
    criteria = submitted.get("criteria")
    if not isinstance(criteria, list):
        return (False, "criteria must be a list (an empty list confirms the deterministic result)")
    rationale = submitted.get("rationale")
    if not isinstance(rationale, str) or not rationale.strip():
        return (False, "rationale must be a non-empty string")
    if PMC_RE.search(rationale):
        return (False, "rationale cites a PMID/literature reference; only evidence-pack values are allowed")
    if not _numbers_allowed(rationale, allowed_numbers):
        return (False, "rationale cites a number not present in the evidence pack")
    det_by_code = {c.code: c for c in det}
    seen: set[str] = set()
    for item in criteria:
        if not isinstance(item, dict):
            return (False, "each criterion must be an object {code, status, strength, reason}")
        code = str(item.get("code", "")).strip()
        if code not in ALL_CRITERION_CODES:
            return (False, f"unknown criterion code {code!r} (must be one of the 28 ACMG codes)")
        if code in seen:
            return (False, f"criterion {code} submitted twice")
        seen.add(code)
        status = str(item.get("status", "")).strip()
        if status not in ALLOWED_STATUSES:
            return (False, f"invalid status {status!r} for {code}")
        strength = item.get("strength")
        if strength is not None and not isinstance(strength, str):
            return (False, f"strength for {code} must be a string or null")
        if strength not in ALLOWED_STRENGTHS:
            return (False, f"invalid strength {strength!r} for {code}")
        base = det_by_code.get(code)
        if base is None:
            return (False, f"criterion {code} is not part of this engine's list")
        changed = (status != base.status) or (strength != base.strength)
        if not changed:
            continue
        if base.status in ("not_evaluable", "not_applicable") and status == "met":
            return (False, f"{code}: a {base.status} criterion may never be upgraded to met")
        if status != "met" and strength is not None:
            return (False, f"{code}: strength applies only to met criteria")
        reason = str(item.get("reason") or "").strip()
        if not reason:
            return (False, f"{code}: any change requires a non-empty reason")
        if PMC_RE.search(reason):
            return (False, f"{code}: reason cites a PMID/literature reference; only evidence-pack values are allowed")
        if not _numbers_allowed(reason, allowed_numbers):
            return (False, f"{code}: reason cites a number not present in the evidence pack")
        if not _mentions_evidence_field(reason):
            return (False, f"{code}: reason must name an evidence field (e.g. {', '.join(EVIDENCE_FIELDS[:8])})")
    return (True, "ok")


def merge_criteria(det: list[Criterion], submitted: dict) -> tuple[list[Criterion], list[str]]:
    """Apply an *already validated* submission to the deterministic list.

    Unchanged criteria keep the engine's entry; changed ones are replaced and
    their codes returned as the adjusted set (reported as adjusted_by_agent).
    """
    by_code = {str(item.get("code")): item for item in (submitted.get("criteria") or [])}
    merged: list[Criterion] = []
    adjusted: list[str] = []
    for c in det:
        item = by_code.get(c.code)
        if item is None:
            merged.append(c)
            continue
        status = str(item.get("status", "")).strip()
        strength = item.get("strength")
        if status == c.status and strength == c.strength:
            merged.append(c)
            continue
        merged.append(Criterion(c.code, status, strength, str(item.get("reason") or "")))
        adjusted.append(c.code)
    return merged, adjusted


# --------------------------------------------------------------------------
# selection and the global run
# --------------------------------------------------------------------------
def _tier_of(cand: dict) -> int:
    try:
        return int(cand.get("tier", 99))
    except (TypeError, ValueError):
        return 99


def _rank_sort_key(cand: dict, idx: int) -> tuple[int, int, int, int]:
    is_lead = 0 if cand.get("is_leading") else 1
    raw_rank = cand.get("rank", idx)
    try:
        r = int(raw_rank)
    except (TypeError, ValueError):
        r = idx
    return (is_lead, r, _tier_of(cand), idx)


def select_review_set(candidates: list[dict], *, max_variants: int = MAX_VARIANTS) -> list[int]:
    """Indices to review: tier 0-2 plus the top 5 overall, capped, in strict rank order (leading first)."""
    order = sorted(range(len(candidates)), key=lambda i: _rank_sort_key(candidates[i], i))
    top5 = set(order[:5])
    eligible = [i for i in order if i in top5 or _tier_of(candidates[i]) <= 2]
    return eligible[:max_variants]


def run_agent(
    ev: Evidence,
    det: list[Criterion],
    backend: Any,
    *,
    max_calls: int = MAX_TOOL_CALLS,
    per_variant_seconds: int = PER_VARIANT_SECONDS,
    max_tokens_note: str = "",
) -> AgentResult:
    """Bounded tool loop for one variant. Never raises: always returns a
    deterministic-safe AgentResult."""
    pack = EvidencePack(ev, det)
    tool_specs = all_tool_specs()
    dispatch = pack.dispatch()
    user_prompt = (
        "Candidate variant evidence pack (JSON):\n"
        f"{build_prompt_context(pack)}\n\n"
        "Review the deterministic criteria above under the contract rules and "
        "submit_assessment. "
        + max_tokens_note
    )
    messages: list[Any] = [Message_system(), Message(role="user", content=user_prompt)]
    deadline = time.monotonic() + per_variant_seconds
    calls = 0

    def fallback(reason: str, status: str) -> AgentResult:
        return AgentResult(criteria=list(det), rationale=reason, status=status,
                           n_tool_calls=calls)

    while True:
        if time.monotonic() >= deadline:
            return fallback("per-variant time budget exhausted before a valid submission", "timeout")
        if calls >= max_calls:
            return fallback("tool-call budget exhausted before a valid submission; deterministic result kept", "rejected")

        try:
            reply = _call_backend_chat(
                backend,
                messages,
                tool_specs,
                trace_meta={"gene": ev.gene, "phase": "review"},
            )
        except Exception as exc:  # backend crash: fail closed, keep engine result
            return fallback(f"LLM backend failed: {exc}", "failed")
        calls += 1

        if not reply.tool_calls:
            messages.append(Message(role="assistant", content=reply.content))
            messages.append(Message(
                role="user",
                content="Call a read tool to inspect the pack, or submit_assessment to finalise.",
            ))
            continue

        messages.append(Message(role="assistant", content=reply.content, tool_calls=list(reply.tool_calls)))
        try:
            accepted = False
            for tc in reply.tool_calls:
                if time.monotonic() >= deadline:
                    return fallback("per-variant time budget exhausted during tool handling", "timeout")
                if calls >= max_calls and tc.name != SUBMIT_TOOL and accepted:
                    continue
                if tc.name == SUBMIT_TOOL:
                    if tc.name not in dispatch or tc.name != SUBMIT_TOOL:
                        pass
                    allowed_numbers = pack_allowed_numbers(pack, det)
                    ok, reason = validate_submission(det, tc.arguments, allowed_numbers=allowed_numbers)
                    if ok:
                        merged, adjusted = merge_criteria(det, tc.arguments)
                        narrative, is_truncated, redacted = _generate_narrative_with_status(backend, ev, merged, deadline)
                        rationale, rat_redacted = strip_recommendations(str(tc.arguments.get("rationale", "")).strip())
                        result = AgentResult(
                            criteria=merged,
                            rationale=rationale or "",
                            status="ok_truncated" if is_truncated else "ok",
                            narrative=narrative,
                            adjusted=adjusted,
                            n_tool_calls=calls,
                            redacted=redacted or rat_redacted,
                        )
                        return result
                    messages.append(Message(role="tool", content=f"SUBMISSION REJECTED: {reason}", name=SUBMIT_TOOL, tool_call_id=tc.id))
                elif tc.name in (set(dispatch) - {SUBMIT_TOOL}):
                    text = dispatch[tc.name](tc.arguments)
                    messages.append(Message(role="tool", content=text, name=tc.name, tool_call_id=tc.id))
                else:
                    messages.append(Message(role="tool", content=f"unknown tool: {tc.name}", name=tc.name, tool_call_id=tc.id))
        except Exception as exc:  # dispatch crash
            return fallback(f"tool dispatch failed: {exc}", "failed")


def Message_system() -> Any:
    return Message(role="system", content=AGENT_SYSTEM_PROMPT)


def pack_allowed_numbers(pack: EvidencePack, det: list[Criterion]) -> set[float]:
    """Numbers legitimately citable: every value in the pack + the engine's
    own rationale strings (thresholds etc. are part of the pack semantics)."""
    allowed: set[float] = set(_extract_numbers(pack.as_json()))
    for c in det:
        allowed |= _extract_numbers(c.rationale)
    return allowed


def _generate_narrative_with_status(
    backend: Any, ev: Evidence, final_criteria: list[Criterion], deadline: float
) -> tuple[str | None, bool, bool]:
    """Narrative plus the recommendation policy.

    Returns ``(narrative, is_truncated, redacted)``. If the narrative recommends a
    clinical action, the model gets one rewrite turn naming the offending
    sentences; if the rewrite still violates (or fails), the sentences are removed.
    """
    text, truncated = _generate_narrative_raw(backend, ev, final_criteria, deadline)
    bad = find_recommendations(text)
    if not bad:
        return text, truncated, False
    if time.monotonic() < deadline:
        try:
            reply = _call_backend_chat(
                backend,
                [
                    Message(role="system", content=NARRATIVE_SYSTEM_PROMPT),
                    Message(role="user", content=f"Narrative for candidate gene {ev.gene}:\n{text}"),
                    Message(role="user", content=violation_message(bad)),
                ],
                tools=[],
                max_tokens=max(int(getattr(backend, "max_tokens", 4096) or 4096) * 2, 8192),
                trace_meta={"gene": ev.gene, "phase": "narrative_policy_retry"},
                attempt="policy_retry",
            )
            text2 = (reply.content or "").strip() if not reply.tool_calls else ""
            if text2 and not find_recommendations(text2) and _ends_with_terminal_punctuation(text2):
                return text2, False, False
        except Exception:
            pass
    clean, _removed = strip_recommendations(text)
    return clean, truncated, True


def _generate_narrative_raw(
    backend: Any, ev: Evidence, final_criteria: list[Criterion], deadline: float
) -> tuple[str | None, bool]:
    """Generate the per-candidate clinical summary paragraph and check for truncation.

    Returns ``(narrative, is_truncated)``. Uses a larger ``max_tokens`` budget
    (at least 8192 tokens, since reasoning tokens count against ``max_tokens``
    on hosted reasoning models), checks both ``finish_reason == 'length'`` and
    terminal punctuation, retries once with an expanded budget when cut, and
    returns ``is_truncated=True`` (which maps to ``agent_status='ok_truncated'``)
    if the retry still fails to produce a complete sentence.
    """
    if time.monotonic() >= deadline:
        return (None, False)
    criteria_table = "\n".join(
        f"- {c.code}: {c.status}" + (f" ({c.strength})" if c.strength else "") + f" — {c.rationale}"
        for c in final_criteria
    )
    user = (
        "Write the report narrative for candidate gene "
        f"{ev.gene}, variant {ev.chrom}:{ev.pos} {ev.ref}>{ev.alt} ({ev.hgvsc or 'hgvs n/a'}). "
        "Base it ONLY on the final criteria list and nothing else:\n"
        f"{criteria_table}"
    )
    base_messages = [
        Message(role="system", content=NARRATIVE_SYSTEM_PROMPT),
        Message(role="user", content=user),
    ]
    base_max = int(getattr(backend, "max_tokens", 4096) or 4096)
    narrative_max_tokens = max(base_max * 2, 8192)
    try:
        reply = _call_backend_chat(
            backend,
            base_messages,
            tools=[],
            max_tokens=narrative_max_tokens,
            trace_meta={"gene": ev.gene, "phase": "narrative"},
            attempt="initial",
        )
    except Exception:
        return (None, False)

    text = (reply.content or "").strip() if not reply.tool_calls else ""
    is_length_cut = getattr(reply, "finish_reason", None) == "length"
    has_term_punct = _ends_with_terminal_punctuation(text)

    if text and not is_length_cut and has_term_punct:
        return (text, False)

    # Retry once if the narrative was truncated (finish_reason=length) or lacks terminal punctuation
    if time.monotonic() < deadline:
        retry_max_tokens = max(narrative_max_tokens * 2, 16384)
        cut_reason = getattr(reply, "finish_reason", None) or "missing_terminal_punctuation"
        if text:
            retry_messages = list(base_messages) + [
                Message(role="assistant", content=text),
                Message(
                    role="user",
                    content=(
                        f"Your previous narrative was cut off before completion ({cut_reason}). "
                        "Write the complete 4-8 sentence clinical summary paragraph from start to finish, "
                        "using only the supplied criteria, and end with a period (.)."
                    ),
                ),
            ]
        else:
            retry_messages = list(base_messages) + [
                Message(
                    role="user",
                    content=(
                        "Write the complete 4-8 sentence clinical summary paragraph from start to finish, "
                        "using only the supplied criteria, and end with a period (.)."
                    ),
                ),
            ]
        try:
            reply2 = _call_backend_chat(
                backend,
                retry_messages,
                tools=[],
                max_tokens=retry_max_tokens,
                trace_meta={"gene": ev.gene, "phase": "narrative_retry"},
                attempt="retry",
            )
            text2 = (reply2.content or "").strip() if not reply2.tool_calls else ""
            is_length_cut2 = getattr(reply2, "finish_reason", None) == "length"
            has_term_punct2 = _ends_with_terminal_punctuation(text2)
            if text2 and not is_length_cut2 and has_term_punct2:
                return (text2, False)
            best = text2 if len(text2) >= len(text) else text
            return (best or None, True)
        except Exception:
            return (text or None, True)

    return (text or None, True)


def _generate_narrative(
    backend: Any, ev: Evidence, final_criteria: list[Criterion], deadline: float
) -> str | None:
    narrative, _truncated, _redacted = _generate_narrative_with_status(backend, ev, final_criteria, deadline)
    return narrative


EXPANSION_SYSTEM_PROMPT = (
    "You are a clinical genetics reviewer writing the grounded interpretation "
    "for a leading short tandem repeat (STR) expansion call from ExpansionHunter "
    "classified with STRchive thresholds.\n"
    "Contract:\n"
    "1. Use get_repeat_call to inspect the RepeatCall fields (gene, alleles/genotype, "
    "CI, STRchive ranges, disease, MOI), then call submit_expansion_narrative with "
    "a concise clinical narrative (3-6 sentences).\n"
    "2. Do NOT adjust or invent ACMG 28-criteria codes; repeat expansions are "
    "evaluated by locus-specific repeat-count thresholds, not SNV/indel ACMG criteria.\n"
    "3. Mention only the repeat expansion gene and values present in the RepeatCall "
    "pack (repeat alleles, confidence interval, normal/intermediate/pathogenic "
    "ranges, mode of inheritance, and disease). Do not cite PMIDs, external "
    "literature, or any number absent from the pack.\n"
    "4. " + RECOMMENDATION_RULE
)


def repeat_pack_allowed_numbers(pack: RepeatPack) -> set[float]:
    """Numbers legitimately citable in a repeat expansion narrative."""
    allowed: set[float] = set(_extract_numbers(pack.as_json()))
    # Allow allele indices 1, 2 and 95 (for '95% CI')
    allowed |= {1.0, 2.0, 95.0}
    return allowed


def validate_expansion_narrative(
    pack: RepeatPack,
    narrative: str,
    *,
    allowed_numbers: set[float] | None = None,
    forbidden_genes: set[str] | None = None,
) -> tuple[bool, str]:
    """Enforce fail-closed grounding contract for a leading repeat expansion narrative."""
    if not isinstance(narrative, str) or not narrative.strip():
        return (False, "narrative must be a non-empty string")
    text = narrative.strip()
    if PMC_RE.search(text):
        return (False, "narrative cites a PMID/literature reference; only RepeatCall pack values are allowed")
    nums = allowed_numbers if allowed_numbers is not None else repeat_pack_allowed_numbers(pack)
    if not _numbers_allowed(text, nums):
        return (False, "narrative cites a number not present in the RepeatCall pack")
    if pack.gene and not re.search(rf"\b{re.escape(pack.gene)}\b", text, re.IGNORECASE):
        return (False, f"narrative must name the repeat expansion gene {pack.gene}")
    for fg in (forbidden_genes or set()):
        if fg and fg.upper() != pack.gene.upper():
            if re.search(rf"\b{re.escape(fg)}\b", text):
                return (False, f"narrative mentions another candidate gene ({fg}) outside the RepeatCall pack")
    return (True, "ok")


def run_expansion_agent(
    repeat_call: Any,
    backend: Any,
    *,
    max_calls: int = MAX_TOOL_CALLS,
    per_variant_seconds: int = PER_VARIANT_SECONDS,
    forbidden_genes: set[str] | None = None,
) -> AgentResult:
    """Bounded, grounded narrative loop for a leading repeat expansion.

    Tool access is restricted to RepeatCall fields only (gene, alleles/genotype,
    CI, STRchive ranges, disease, MOI). No ACMG criteria adjustments are allowed.
    """
    pack = RepeatPack(repeat_call)
    allowed_numbers = repeat_pack_allowed_numbers(pack)
    tool_specs = repeat_tool_specs()
    user_prompt = (
        "Leading repeat expansion RepeatCall pack (JSON):\n"
        f"{pack.as_json()}\n\n"
        "Inspect the RepeatCall fields (gene, alleles, CI, STRchive ranges, disease, MOI) "
        "using get_repeat_call if needed, and submit a grounded clinical narrative via "
        "submit_expansion_narrative. Do not adjust any ACMG criteria."
    )
    messages: list[Any] = [
        Message(role="system", content=EXPANSION_SYSTEM_PROMPT),
        Message(role="user", content=user_prompt),
    ]
    deadline = time.monotonic() + per_variant_seconds
    calls = 0

    policy_retry_used = False

    def fallback(reason: str, status: str) -> AgentResult:
        return AgentResult(
            criteria=[],
            rationale=reason,
            status=status,
            narrative=None,
            adjusted=[],
            n_tool_calls=calls,
        )

    while True:
        if time.monotonic() >= deadline:
            return fallback("per-variant time budget exhausted before a valid expansion narrative", "timeout")
        if calls >= max_calls:
            return fallback("tool-call budget exhausted before a valid expansion narrative; deterministic result kept", "rejected")

        try:
            reply = _call_backend_chat(
                backend,
                messages,
                tool_specs,
                trace_meta={"gene": pack.gene, "phase": "expansion"},
            )
        except Exception as exc:
            return fallback(f"LLM backend failed: {exc}", "failed")
        calls += 1

        if not reply.tool_calls:
            content_text = (reply.content or "").strip()
            if content_text:
                ok, reason = validate_expansion_narrative(
                    pack,
                    content_text,
                    allowed_numbers=allowed_numbers,
                    forbidden_genes=forbidden_genes,
                )
                bad = find_recommendations(content_text) if ok else []
                if ok and bad and not policy_retry_used:
                    policy_retry_used = True
                    ok, reason = False, violation_message(bad)
                if ok:
                    clean, red = strip_recommendations(content_text)
                    return AgentResult(
                        criteria=[],
                        rationale=clean or "",
                        status="ok",
                        narrative=clean,
                        adjusted=[],
                        n_tool_calls=calls,
                        redacted=red,
                    )
                messages.append(Message(role="assistant", content=reply.content))
                messages.append(Message(
                    role="user",
                    content=f"NARRATIVE REJECTED: {reason}. Call get_repeat_call or submit_expansion_narrative with a grounded narrative for {pack.gene}.",
                ))
            else:
                messages.append(Message(role="assistant", content=reply.content))
                messages.append(Message(
                    role="user",
                    content=f"Call get_repeat_call to inspect the pack, or submit_expansion_narrative to finalise the narrative for {pack.gene}.",
                ))
            continue

        messages.append(Message(role="assistant", content=reply.content, tool_calls=list(reply.tool_calls)))
        try:
            for tc in reply.tool_calls:
                if time.monotonic() >= deadline:
                    return fallback("per-variant time budget exhausted during tool handling", "timeout")
                if tc.name in (REPEAT_SUBMIT_TOOL, SUBMIT_TOOL):
                    args_dict = tc.arguments if isinstance(tc.arguments, dict) else {}
                    if args_dict.get("criteria"):
                        messages.append(Message(
                            role="tool",
                            content="SUBMISSION REJECTED: repeat expansion review does not allow ACMG criteria adjustments",
                            name=tc.name,
                            tool_call_id=tc.id,
                        ))
                        continue
                    nar = str(args_dict.get("narrative") or args_dict.get("rationale") or "").strip()
                    ok, reason = validate_expansion_narrative(
                        pack,
                        nar,
                        allowed_numbers=allowed_numbers,
                        forbidden_genes=forbidden_genes,
                    )
                    bad = find_recommendations(nar) if ok else []
                    if ok and bad and not policy_retry_used:
                        policy_retry_used = True
                        ok, reason = False, violation_message(bad)
                    if ok:
                        clean, red = strip_recommendations(nar)
                        rat, red2 = strip_recommendations(str(args_dict.get("rationale") or nar).strip())
                        return AgentResult(
                            criteria=[],
                            rationale=rat or "",
                            status="ok",
                            narrative=clean,
                            adjusted=[],
                            n_tool_calls=calls,
                            redacted=red or red2,
                        )
                    messages.append(Message(
                        role="tool",
                        content=f"SUBMISSION REJECTED: {reason}",
                        name=tc.name,
                        tool_call_id=tc.id,
                    ))
                elif tc.name == REPEAT_READ_TOOL:
                    messages.append(Message(
                        role="tool",
                        content=pack.as_json(),
                        name=tc.name,
                        tool_call_id=tc.id,
                    ))
                else:
                    messages.append(Message(
                        role="tool",
                        content=f"unknown tool for repeat expansion: {tc.name}",
                        name=tc.name,
                        tool_call_id=tc.id,
                    ))
        except Exception as exc:
            return fallback(f"tool dispatch failed: {exc}", "failed")


def run_all(
    candidates: list[dict],
    backend: Any,
    *,
    global_seconds: int = GLOBAL_SECONDS,
    max_variants: int = MAX_VARIANTS,
    per_variant_seconds: int = PER_VARIANT_SECONDS,
    start_time: float | None = None,
    leading_expansion: dict | None = None,
) -> tuple[list[dict], str]:
    """Review the selected candidates in-place; return (candidates, global status).

    Selection and budgets: tier 0-2 + top 5 overall, ≤ max_variants, per-variant
    12 tool calls / per_variant_seconds, global_seconds wall clock. Anyone not
    reviewed keeps the deterministic criteria with agent_status="skipped" and
    the caller reports llm.status ("ok"/"partial"/...) from the return value.

    When ``leading_expansion`` is provided (e.g. pathogenic-range HTT expansion),
    it is reviewed FIRST with the full per_variant_seconds budget via
    ``run_expansion_agent`` before any SNV/indel candidate consumes the global
    budget.

    Each candidate dict is mutated: criteria -> final list, plus keys
    agent_status, agent_rationale, narrative, adjusted.
    """
    started = start_time if start_time is not None else time.monotonic()
    deadline = started + global_seconds
    selected = select_review_set(candidates, max_variants=max_variants)
    for cand in candidates:
        cand["agent_status"] = "skipped"
        cand["agent_rationale"] = None
        cand["narrative"] = None
        cand["adjusted"] = []

    if leading_expansion is not None:
        leading_expansion["agent_status"] = "skipped"
        leading_expansion["agent_rationale"] = None
        leading_expansion["narrative"] = None
        leading_expansion["agent_narrative"] = None
        leading_expansion["adjusted"] = []
        leading_expansion["adjusted_by_agent"] = []

    if backend is None:
        return candidates, OFF_STATUS
    if not selected and leading_expansion is None:
        return candidates, GLOBAL_STATUS_OK

    reviewed = fallbacks = 0
    exp_ok = True
    if leading_expansion is not None:
        t_exp = time.monotonic()
        exp_budget = int(per_variant_seconds) if per_variant_seconds <= 0 else max(1, int(per_variant_seconds))
        snv_genes = {
            str(c.get("gene") or getattr(c.get("evidence"), "gene", ""))
            for c in candidates
        } - {str(leading_expansion.get("gene") or "")}
        exp_res = run_expansion_agent(
            leading_expansion,
            backend,
            per_variant_seconds=exp_budget,
            forbidden_genes=snv_genes,
        )
        leading_expansion["wall_seconds"] = round(time.monotonic() - t_exp, 2)
        leading_expansion["n_tool_calls"] = exp_res.n_tool_calls
        leading_expansion["agent_status"] = exp_res.status
        leading_expansion["agent_rationale"] = exp_res.rationale
        leading_expansion["narrative"] = exp_res.narrative if exp_res.status in ("ok", "ok_truncated") else None
        leading_expansion["agent_narrative"] = exp_res.narrative if exp_res.status in ("ok", "ok_truncated") else None
        leading_expansion["adjusted"] = []
        leading_expansion["adjusted_by_agent"] = []
        leading_expansion["narrative_redacted"] = bool(exp_res.redacted)
        exp_ok = exp_res.status == "ok"
        if not exp_ok:
            fallbacks += 1

    lead_idx = selected[0] if selected else None
    for pos, i in enumerate(selected):
        cand = candidates[i]
        is_first_overall = (pos == 0 and leading_expansion is None)
        # The leading finding always runs first and gets the full
        # per_variant_seconds budget before any subsequent candidate can consume
        # the global budget.
        if not is_first_overall and time.monotonic() >= deadline:
            cand["agent_status"] = "skipped"
            fallbacks += 1
            continue
        reviewed += 1
        t_var = time.monotonic()
        if is_first_overall:
            var_budget = int(per_variant_seconds) if per_variant_seconds <= 0 else max(1, int(per_variant_seconds))
        else:
            if per_variant_seconds <= 0:
                var_budget = 0
            else:
                remaining_global = deadline - t_var
                var_budget = max(1, int(min(per_variant_seconds, max(1, remaining_global))))
        result = run_agent(
            cand["evidence"],
            cand["criteria"],
            backend,
            per_variant_seconds=var_budget,
        )
        cand["wall_seconds"] = round(time.monotonic() - t_var, 2)
        cand["n_tool_calls"] = result.n_tool_calls
        cand["agent_status"] = result.status
        cand["agent_rationale"] = result.rationale
        cand["narrative"] = result.narrative
        cand["adjusted"] = list(result.adjusted)
        cand["criteria"] = result.criteria
        cand["narrative_redacted"] = bool(result.redacted)
        if result.status != "ok":
            fallbacks += 1

    if leading_expansion is not None:
        lead_ok = exp_ok
    else:
        lead_ok = (lead_idx is not None and candidates[lead_idx].get("agent_status") == "ok")

    if not lead_ok:
        return candidates, GLOBAL_STATUS_PARTIAL
    if fallbacks == 0 and reviewed == len(selected):
        return candidates, GLOBAL_STATUS_OK
    return candidates, GLOBAL_STATUS_PARTIAL


def run_llm_layer(
    candidates: list[dict],
    backend: Any,
    *,
    global_seconds: int = GLOBAL_SECONDS,
    max_variants: int = MAX_VARIANTS,
    per_variant_seconds: int = PER_VARIANT_SECONDS,
) -> list[dict]:
    """Scope-note entry point: same dicts, criteria replaced by the adjusted list
    plus agent_rationale / agent_status / narrative. Use run_all() when the
    caller also needs the run-level llm.status."""
    result, _status = run_all(
        candidates, backend,
        global_seconds=global_seconds,
        max_variants=max_variants,
        per_variant_seconds=per_variant_seconds,
    )
    return result
