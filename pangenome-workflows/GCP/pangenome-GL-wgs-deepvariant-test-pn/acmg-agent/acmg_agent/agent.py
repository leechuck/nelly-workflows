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

from acmg_agent.llm_backends import Message  # noqa: F401 (re-export for callers)
from acmg_agent.model import Criterion, Evidence
from acmg_agent.tools import (
    SUBMIT_TOOL,
    EvidencePack,
    all_tool_specs,
    build_prompt_context,
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
    "If you cannot improve the engine's assessment with pack evidence, submit "
    "no changes and say so in the rationale."
)


@dataclass
class AgentResult:
    """Outcome of one per-variant review loop."""
    criteria: list[Criterion]
    rationale: str = ""
    status: str = "ok"            # ok | rejected | failed | timeout
    narrative: str | None = None
    adjusted: list[str] = field(default_factory=list)
    n_tool_calls: int = 0


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
    "in neutral terms."
)


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


def select_review_set(candidates: list[dict], *, max_variants: int = MAX_VARIANTS) -> list[int]:
    """Indices to review: tier 0-2 plus the top 5 overall, capped, rank order."""
    order = sorted(range(len(candidates)), key=lambda i: candidates[i].get("rank", 10 ** 6))
    top5 = set(order[:5])
    tiers = {i for i in order if _tier_of(candidates[i]) <= 2}
    selected = set(sorted(tiers | top5)[:max_variants])
    return [i for i in order if i in selected]


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
            reply = backend.chat(messages, tool_specs)
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
                        result = AgentResult(
                            criteria=merged,
                            rationale=str(tc.arguments.get("rationale", "")).strip(),
                            status="ok",
                            adjusted=adjusted,
                            n_tool_calls=calls,
                        )
                        result.narrative = _generate_narrative(backend, ev, merged, deadline)
                        return result
                    messages.append(Message(role="tool", content=f"SUBMISSION REJECTED: {reason}", name=SUBMIT_TOOL))
                elif tc.name in (set(dispatch) - {SUBMIT_TOOL}):
                    text = dispatch[tc.name](tc.arguments)
                    messages.append(Message(role="tool", content=text, name=tc.name))
                else:
                    messages.append(Message(role="tool", content=f"unknown tool: {tc.name}", name=tc.name))
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


def _generate_narrative(
    backend: Any, ev: Evidence, final_criteria: list[Criterion], deadline: float
) -> str | None:
    if time.monotonic() >= deadline:
        return None
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
    try:
        reply = backend.chat(
            [Message(role="system", content=NARRATIVE_SYSTEM_PROMPT), Message(role="user", content=user)],
            tools=[],
        )
    except Exception:
        return None
    if not reply.tool_calls and reply.content and reply.content.strip():
        return reply.content.strip()
    return None


def run_all(
    candidates: list[dict],
    backend: Any,
    *,
    global_seconds: int = GLOBAL_SECONDS,
    max_variants: int = MAX_VARIANTS,
    per_variant_seconds: int = PER_VARIANT_SECONDS,
    start_time: float | None = None,
) -> tuple[list[dict], str]:
    """Review the selected candidates in-place; return (candidates, global status).

    Selection and budgets: tier 0-2 + top 5 overall, ≤ max_variants, per-variant
    12 tool calls / per_variant_seconds, global_seconds wall clock. Anyone not
    reviewed keeps the deterministic criteria with agent_status="skipped" and
    the caller reports llm.status ("ok"/"partial"/...) from the return value.

    Each candidate dict is mutated: criteria -> final list, plus keys
    agent_status, agent_rationale, narrative, adjusted.
    """
    started = start_time if start_time is not None else time.monotonic()
    deadline = started + global_seconds
    selected = select_review_set(candidates, max_variants=max_variants)
    reviewed = fallbacks = 0
    for i, cand in enumerate(candidates):
        if i not in selected:
            cand["agent_status"] = "skipped"
            cand["agent_rationale"] = None
            cand["narrative"] = None
            cand["adjusted"] = []
            continue
        if backend is None:
            cand["agent_status"] = "skipped"
            cand["agent_rationale"] = None
            cand["narrative"] = None
            cand["adjusted"] = []
            continue
        if time.monotonic() >= deadline:
            cand["agent_status"] = "skipped"
            cand["agent_rationale"] = None
            cand["narrative"] = None
            cand["adjusted"] = []
            continue
        reviewed += 1
        remaining_global = deadline - time.monotonic()
        budget = min(per_variant_seconds, max(1, remaining_global))
        result = run_agent(
            cand["evidence"], cand["criteria"], backend,
            per_variant_seconds=max(1, int(budget)),
        )
        cand["agent_status"] = result.status
        cand["agent_rationale"] = result.rationale
        cand["narrative"] = result.narrative
        cand["adjusted"] = list(result.adjusted)
        cand["criteria"] = result.criteria
        if result.status != "ok":
            fallbacks += 1
    if not selected or reviewed == 0:
        return candidates, OFF_STATUS if backend is None else GLOBAL_STATUS_OK
    if fallbacks == 0:
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
