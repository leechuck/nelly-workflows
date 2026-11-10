"""Tests for the bounded LLM agent layer (FakeBackend only -- no real model)."""
from __future__ import annotations

from typing import Any

from acmg_agent.agent import (
    run_agent,
    run_all,
    run_expansion_agent,
    run_llm_layer,
    select_review_set,
    validate_expansion_narrative,
    validate_submission,
)
from acmg_agent.criteria import evaluate_all
from acmg_agent.llm_backends import (
    Message,
    ToolCall,
    select_backend_from_env,
)
from acmg_agent.model import Criterion, Evidence


def make_ev(**overrides) -> Evidence:
    base: dict[str, Any] = {
        "gene": "TEST", "moi": "AR", "chrom": "1", "pos": 100, "ref": "A", "alt": "G",
        "zygosity": "hom", "transcript": "NM_000001.1", "mane_select": True,
        "consequence": "splice_donor_variant", "hgvsc": "c.100+1G>A", "hgvsp": "p.?",
        "exon": "2/10", "intron": "1/9", "cds_position": 100, "protein_position": None,
        "nmd_escape": False, "in_repeat_region": False,
        "max_af": 0.0, "popmax_af": 0.0, "nhomalt": 0, "bs1_threshold": 0.005,
        "revel": 0.5, "spliceai": 0.05,
        "clinvar_id": None, "clinvar_sig": None, "clinvar_stars": 0,
        "same_aa_ClinVar_PLP": 0, "diff_aa_ClinVar_PLP": 0,
        "lof_mechanism": True, "loeuf": 0.2, "pli": 0.95,
        "clinvar_missense_plp": 5, "clinvar_missense_blb": 10, "in_hotspot": False,
        "other_var_in_gene": [], "human_pheno_score": 0.7, "pheno_evidence": "HPO:0001",
        "disease": "TEST disease",
    }
    base.update(overrides)
    return Evidence(**base)


def det_criteria(ev: Evidence) -> list[Criterion]:
    return evaluate_all(ev)


class FakeBackend:
    """Scripted tool-calling backend.

    Script items: ToolCall (one call), Message (custom message with finish_reason/usage),
    str (plain assistant content), Exception instance (raised). When the script is
    exhausted, a canned content-only reply ending with a period is returned.
    """

    def __init__(self, script):
        self.script = list(script)
        self.calls: list[tuple[list[Message], list[dict]]] = []
        self.max_tokens_seen: list[int | None] = []
        self.trace_path: str | None = None

    def chat(
        self,
        messages: list[Message],
        tools: list[dict],
        *,
        max_tokens: int | None = None,
        trace_meta: dict | None = None,
    ) -> Message:
        self.calls.append((messages, tools))
        self.max_tokens_seen.append(max_tokens)
        if not self.script:
            return Message(role="assistant", content=f"narrative for candidate {len(self.calls)}.", finish_reason="stop")
        item = self.script.pop(0)
        if isinstance(item, Exception):
            raise item
        if isinstance(item, Message):
            return item
        if isinstance(item, str):
            return Message(role="assistant", content=item, finish_reason="stop")
        return Message(role="assistant", tool_calls=[item], finish_reason="tool_calls")


def submit_tool(payload: dict) -> ToolCall:
    return ToolCall(name="submit_assessment", arguments=payload)


# ---------------------------------------------------------------- Run agent --
def test_agent_accepts_strength_downgrade():
    ev = make_ev()
    det = det_criteria(ev)
    backend = FakeBackend([
        ToolCall(name="get_variant_evidence", arguments={}),
        submit_tool({
            "criteria": [
                {"code": "PVS1", "status": "met", "strength": "strong",
                 "reason": "nmd_escape is True in VEP, PVS1_Strong per pack"},
            ],
            "rationale": "PVS1 is met, downgraded to strong (reason in criteria).",
        }),
    ])
    result = run_agent(ev, det, backend)
    assert result.status == "ok"
    pvs1 = [c for c in result.criteria if c.code == "PVS1"][0]
    assert pvs1.status == "met" and pvs1.strength == "strong"
    assert "nmd_escape" in pvs1.rationale
    assert result.adjusted == ["PVS1"]
    # untouched criteria keep the deterministic entry
    ps2_final = [c for c in result.criteria if c.code == "PS2"][0]
    ps2_det = [c for c in det if c.code == "PS2"][0]
    assert ps2_final == ps2_det


def test_agent_rejects_not_evaluable_to_met():
    ev = make_ev()
    det = det_criteria(ev)
    assert [c for c in det if c.code == "PS2"][0].status == "not_evaluable"
    backend = FakeBackend([
        submit_tool({
            "criteria": [
                {"code": "PS2", "status": "met", "strength": "strong",
                 "reason": "zygosity is hom so the variant must be de novo"},
            ],
            "rationale": "trying to upgrade a not_evaluable criterion",
        }),
        submit_tool({
            "criteria": [
                {"code": "PVS1", "status": "met", "strength": "moderate",
                 "reason": "exon is 2/10 per pack, downgraded"},
            ],
            "rationale": "safe retry submission",
        }),
    ])
    result = run_agent(ev, det, backend)
    assert result.status == "ok"
    assert [c for c in result.criteria if c.code == "PS2"][0].status == "not_evaluable"
    # the PS2 met attempt must never survive
    assert all(c.code != "PS2" or c.status != "met" for c in result.criteria)


def test_agent_rejects_unknown_code():
    ev = make_ev()
    det = det_criteria(ev)
    backend = FakeBackend([
        submit_tool({
            "criteria": [
                {"code": "PVS99", "status": "met", "strength": "very_strong",
                 "reason": "invented code"},
            ],
            "rationale": "inventing a code",
        }),
    ])
    result = run_agent(ev, det, backend)
    assert result.status == "rejected"
    assert [c for c in result.criteria if c.code == "PVS99"] == []
    assert result.criteria == det            # full deterministic fallback
    assert result.narrative is None

def test_agent_budget_exhausted_falls_back():
    ev = make_ev()
    det = det_criteria(ev)
    # model that keeps pinging read tools, never submits
    script = [ToolCall(name="get_gene_summary", arguments={})] * 20
    backend = FakeBackend(script)
    result = run_agent(ev, det, backend, max_calls=3)
    assert result.status == "rejected"
    assert result.n_tool_calls == 3
    assert result.criteria == det


def test_agent_backend_exception_status_failed():
    ev = make_ev()
    det = det_criteria(ev)
    result = run_agent(ev, det, FakeBackend([RuntimeError("llama exploded")]))
    assert result.status == "failed"
    assert result.criteria == det


def test_agent_timeout_status():
    ev = make_ev()
    det = det_criteria(ev)
    result = run_agent(ev, det, FakeBackend([]), per_variant_seconds=0)
    assert result.status == "timeout"
    assert result.criteria == det


def test_agent_pmid_or_pack_foreign_number_rejected():
    ev = make_ev()
    det = det_criteria(ev)
    # threshold 0.05 is not a value in this pack (max_af=0.0 ...)
    backend = FakeBackend([
        submit_tool({
            "criteria": [
                {"code": "PVS1", "status": "met", "strength": "strong",
                 "reason": "Abou Tayoun 2018 PMID 25741868 says use 0.05"},
            ],
            "rationale": "citing literature",
        }),
    ])
    result = run_agent(ev, det, backend)
    assert result.status == "rejected"
    assert result.criteria == det


# -------------------------------------------------------- validate_submission --
def test_validator_change_requires_evidence_field_name():
    ev = make_ev()
    det = det_criteria(ev)
    ok, _ = validate_submission(det, {
        "criteria": [{"code": "PVS1", "status": "met", "strength": "moderate",
                      "reason": "looks weaker to me"}],
        "rationale": "downgrade",
    })
    assert not ok
    ok, reason = validate_submission(det, {
        "criteria": [{"code": "PVS1", "status": "met", "strength": "moderate",
                      "reason": "nmd_escape in the pack points at a leader exon"}],
        "rationale": "downgrade",
    })
    assert ok, reason


def test_validator_unchanged_criteria_may_be_omitted():
    ev = make_ev()
    det = det_criteria(ev)
    ok, reason = validate_submission(det, {
        "criteria": [{"code": "PVS1", "status": "met", "strength": "very_strong", "reason": ""}],
        "rationale": "no changes, confirming the deterministic result",
    })
    assert ok, reason


def test_validator_strength_only_on_met():
    ev = make_ev()
    det = det_criteria(ev)
    ok, _ = validate_submission(det, {
        "criteria": [{"code": "PM2", "status": "not_met", "strength": "supporting",
                      "reason": "max_af is not None here"}],
        "rationale": "bad shape",
    })
    assert not ok


# ------------------------------------------------------- multi-candidate run --
def _candidate(rank: int, tier: int, gene: str) -> dict:
    ev = make_ev(gene=gene)
    return {"rank": rank, "tier": tier, "gene": gene, "evidence": ev,
            "criteria": det_criteria(ev)}


def test_run_llm_layer_selection_tier02_top5_cap8():
    cands = [_candidate(r, 3 if r > 3 else r - 1, f"GENE{r}") for r in range(1, 6)]
    cands += [_candidate(r, 3, f"GENE{r}") for r in range(6, 11)]
    backend = FakeBackend([])   # canned content replies -> immediate rejected fallbacks
    out = run_llm_layer(cands, backend, per_variant_seconds=0)
    statuses = {c["gene"]: c["agent_status"] for c in out}
    # ranks 1-3 carry tiers 0-2, ranks 4-5 are in the top5 overall: reviewed.
    # GENE6..GENE10 are not tier 0-2 and not in the top 5: skipped.
    for gene in ("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"):
        assert statuses[gene] != "skipped"
    for gene in ("GENE6", "GENE7", "GENE8", "GENE9", "GENE10"):
        assert statuses[gene] == "skipped"


def test_run_llm_layer_capEight():
    cands = [_candidate(r, 3, f"G{r}") for r in range(1, 10)]
    assert len(select_review_set(cands)) == 5


def test_run_llm_layer_global_partial_when_budget_left_out():
    cands = [_candidate(r, r - 1 if r <= 3 else 3, f"G{r}") for r in range(1, 10)]
    backend = FakeBackend([])
    _out, status = run_all(cands, backend, per_variant_seconds=0, global_seconds=60)
    assert status == "partial"


def test_run_llm_layer_global_ok_when_all_reviewed_ok():
    cands = [_candidate(1, 0, "AAA")]
    payload = {
        "criteria": [
            {"code": "PVS1", "status": "met", "strength": "strong",
             "reason": "nmd_escape from the pack"},
        ],
        "rationale": "downgraded one strength",
    }
    backend = FakeBackend([submit_tool(payload)])
    out, status = run_all(cands, backend)
    assert status == "ok"
    assert out[0]["agent_status"] == "ok"
    assert out[0]["adjusted"] == ["PVS1"]
    pvs1 = [c for c in out[0]["criteria"] if c.code == "PVS1"][0]
    assert pvs1.strength == "strong"


def test_run_llm_layer_backend_none_all_skipped_off():
    cands = [_candidate(1, 0, "AAA")]
    out, status = run_all(cands, None)
    assert status == "off"
    assert out[0]["agent_status"] == "skipped"
    assert out[0]["criteria"] == det_criteria(out[0]["evidence"])


def test_narrative_mentions_only_candidate_genes():
    cands = [_candidate(1, 0, "HEXA"), _candidate(2, 3, "HEXB")]
    payload = {
        "criteria": [
            {"code": "PVS1", "status": "met", "strength": "strong",
             "reason": "nmd_escape from the pack"},
        ],
        "rationale": "HEXA downgrade per nmd_escape",
    }
    backend = FakeBackend([
        submit_tool(payload),
        "HEXA is described as casually affected by the splice consequence; "
        "a nearby near-tie remains in the differential.",
    ])
    out, status = run_all(cands, backend)
    # HEXB has no scripted submission so it falls back (run -> partial); the
    # HEXA review itself succeeded.
    assert out[0]["agent_status"] == "ok"
    narrative = out[0]["narrative"]
    assert narrative and "HEXA" in narrative
    # the narrative names no gene outside the candidate set
    assert "HEXB" not in narrative


def test_run_all_reviews_in_rank_order_starting_with_leading():
    import re
    # Input list is out of rank order: rank 3 first, leading (rank 1) in middle, rank 2 last
    c_rank3 = _candidate(3, 2, "GENE3")
    c_lead = _candidate(1, 0, "HEXA")
    c_lead["is_leading"] = True
    c_rank2 = _candidate(2, 1, "SEMA6B")
    cands = [c_rank3, c_lead, c_rank2]

    payload = {
        "criteria": [
            {"code": "PVS1", "status": "met", "strength": "strong",
             "reason": "nmd_escape from the pack"},
        ],
        "rationale": "ok review",
    }
    backend = FakeBackend([
        submit_tool(payload), "HEXA narrative.",
        submit_tool(payload), "SEMA6B narrative.",
        submit_tool(payload), "GENE3 narrative.",
    ])
    out, status = run_all(cands, backend)
    assert status == "ok"
    # Extract the order in which genes were first presented to backend.chat
    seen_genes: list[str] = []
    for msgs, _tools in backend.calls:
        user_txt = next((m.content or "" for m in msgs if m.role == "user"), "")
        m = re.search(r'"gene":\s*"([^"]+)"', user_txt)
        if m and (not seen_genes or seen_genes[-1] != m.group(1)):
            seen_genes.append(m.group(1))
    assert seen_genes == ["HEXA", "SEMA6B", "GENE3"]


def test_run_all_gives_leading_full_budget_even_when_global_tight():
    import time
    c_lead = _candidate(1, 0, "HEXA")
    c_lead["is_leading"] = True
    c_other = _candidate(2, 1, "SEMA6B")
    cands = [c_other, c_lead]

    payload = {
        "criteria": [
            {"code": "PVS1", "status": "met", "strength": "strong",
             "reason": "nmd_escape from the pack"},
        ],
        "rationale": "leading reviewed despite global elapsed",
    }
    backend = FakeBackend([submit_tool(payload), "HEXA leading narrative."])
    # Global budget already exhausted before run_all starts, but leading candidate
    # must still receive its full per_variant_seconds budget first.
    out, status = run_all(
        cands,
        backend,
        per_variant_seconds=60,
        global_seconds=10,
        start_time=time.monotonic() - 100,
    )
    by_gene = {c["gene"]: c for c in out}
    assert by_gene["HEXA"]["agent_status"] == "ok"
    assert by_gene["SEMA6B"]["agent_status"] == "skipped"
    assert status == "partial"


def test_run_all_never_ok_when_leading_candidate_unreviewed():
    c_lead = _candidate(1, 0, "HEXA")
    c_lead["is_leading"] = True
    c_other = _candidate(2, 1, "SEMA6B")
    cands = [c_lead, c_other]

    ok_payload = {
        "criteria": [
            {"code": "PVS1", "status": "met", "strength": "strong",
             "reason": "nmd_escape from the pack"},
        ],
        "rationale": "lower candidate ok",
    }
    # First candidate (HEXA) fails with RuntimeError; second candidate (SEMA6B) succeeds
    backend = FakeBackend([
        RuntimeError("timeout on leading"),
        submit_tool(ok_payload),
        "SEMA6B narrative.",
    ])
    out, status = run_all(cands, backend)
    by_gene = {c["gene"]: c for c in out}
    assert by_gene["HEXA"]["agent_status"] == "failed"
    assert by_gene["SEMA6B"]["agent_status"] == "ok"
    assert status == "partial"


def _make_htt_expansion() -> dict:
    return {
        "gene": "HTT",
        "genotypes": "18/46",
        "ci": "18-18/45-47",
        "allele1": 18,
        "allele2": 46,
        "max_allele": 46,
        "normal_max": 26,
        "pathogenic_min": 36,
        "moi": "AD",
        "disease": "Huntington disease (OMIM:143100)",
        "status": "pathogenic_range",
        "tier": "pathogenic",
        "ranges": {"benign": [1, 26], "intermediate": [27, 35], "pathogenic": [36, 250]},
    }


def test_expansion_agent_accepts_grounded_narrative_with_repeat_call_tool():
    exp = _make_htt_expansion()
    narrative_txt = (
        "HTT carries a pathogenic-range repeat expansion with alleles 18/46 "
        "(95% CI 18-18/45-47), where the expanded allele of 46 repeats exceeds "
        "the pathogenic threshold of 36 repeats (normal <= 26) for autosomal "
        "dominant (AD) Huntington disease (OMIM:143100)."
    )
    backend = FakeBackend([
        ToolCall(name="get_repeat_call", arguments={}),
        ToolCall(name="submit_expansion_narrative", arguments={"narrative": narrative_txt}),
    ])
    res = run_expansion_agent(exp, backend, forbidden_genes={"ATXN1", "CACNA1A"})
    assert res.status == "ok"
    assert res.adjusted == []
    assert res.criteria == []
    assert res.narrative == narrative_txt
    assert res.n_tool_calls == 2


def test_expansion_agent_rejects_criteria_adjustments_and_ungrounded_numbers():
    exp = _make_htt_expansion()
    backend = FakeBackend([
        # Attempt 1: tries to adjust ACMG criteria -> rejected
        ToolCall(
            name="submit_assessment",
            arguments={
                "criteria": [{"code": "PVS1", "status": "met", "strength": "very_strong", "reason": "repeat"}],
                "rationale": "HTT expansion",
            },
        ),
        # Attempt 2: cites ungrounded number 999 -> rejected
        ToolCall(
            name="submit_expansion_narrative",
            arguments={"narrative": "HTT expansion of 999 repeats causes disease."},
        ),
        # Attempt 3: cites PMID -> rejected
        ToolCall(
            name="submit_expansion_narrative",
            arguments={"narrative": "HTT expansion of 46 repeats per PMID 12345678."},
        ),
        # Attempt 4: names forbidden gene ATXN1 -> rejected
        ToolCall(
            name="submit_expansion_narrative",
            arguments={"narrative": "HTT and ATXN1 have 46 repeats."},
        ),
    ])
    res = run_expansion_agent(exp, backend, max_calls=4, forbidden_genes={"ATXN1"})
    assert res.status == "rejected"
    assert res.narrative is None


def test_run_all_reviews_leading_expansion_first_and_requires_it_for_ok():
    import time
    exp = _make_htt_expansion()
    c_snv = _candidate(1, 2, "ATXN1")
    narrative_txt = "HTT repeat expansion 18/46 (95% CI 18-18/45-47) >= 36 repeats for AD Huntington disease."
    backend = FakeBackend([
        ToolCall(name="get_repeat_call", arguments={}),
        ToolCall(name="submit_expansion_narrative", arguments={"narrative": narrative_txt}),
    ])
    # Even when global budget is already elapsed, leading_expansion gets full per_variant_seconds first
    out, status = run_all(
        [c_snv],
        backend,
        per_variant_seconds=60,
        global_seconds=10,
        start_time=time.monotonic() - 100,
        leading_expansion=exp,
    )
    assert exp["agent_status"] == "ok"
    assert exp["agent_narrative"] == narrative_txt
    assert out[0]["agent_status"] == "skipped"
    assert status == "partial"

    # When leading_expansion fails even though SNV succeeds, global status is still partial
    exp2 = _make_htt_expansion()
    ok_snv_payload = {
        "criteria": [{"code": "PVS1", "status": "met", "strength": "strong", "reason": "nmd_escape in pack"}],
        "rationale": "ok snv",
    }
    backend2 = FakeBackend([
        RuntimeError("expansion timeout"),
        submit_tool(ok_snv_payload),
        "ATXN1 narrative.",
    ])
    out2, status2 = run_all([_candidate(1, 2, "ATXN1")], backend2, leading_expansion=exp2)
    assert exp2["agent_status"] == "failed"
    assert out2[0]["agent_status"] == "ok"
    assert status2 == "partial"


# ----------------------------------------------------------------- backends --
def test_select_backend_from_env(monkeypatch):
    monkeypatch.setenv("ACMG_LLM_BACKEND", "off")
    assert select_backend_from_env() is None
    monkeypatch.setenv("ACMG_LLM_BACKEND", "llama")
    monkeypatch.delenv("ACMG_MODEL_PATH", raising=False)
    assert select_backend_from_env() is None
    monkeypatch.setenv("ACMG_LLM_BACKEND", "")
    assert select_backend_from_env() is None


def test_message_conversion_roundtrip():
    from acmg_agent.llm_backends import _completion_to_message, _message_to_openai
    m = Message(role="assistant", tool_calls=[ToolCall(name="t", arguments={"a": 1})])
    d = _message_to_openai(m)
    assert d["tool_calls"][0]["function"]["arguments"] == '{"a": 1}'
    m2 = _completion_to_message({"choices": [{"message": d}]})
    assert m2.tool_calls[0].name == "t" and m2.tool_calls[0].arguments == {"a": 1}


def test_openai_backend_max_tokens_and_reasoning_ignored(monkeypatch):
    from acmg_agent.llm_backends import OpenAICompatBackend, _completion_to_message, _message_to_openai

    monkeypatch.delenv("ACMG_LLM_MAX_TOKENS", raising=False)
    b1 = OpenAICompatBackend("http://unimatrix01.kaust.edu.sa:8000/v1", "qwen3.8-27b")
    assert b1.max_tokens == 4096
    assert b1.timeout_seconds == 300

    monkeypatch.setenv("ACMG_LLM_MAX_TOKENS", "8192")
    b2 = OpenAICompatBackend("http://unimatrix01.kaust.edu.sa:8000/v1", "qwen3.8-27b")
    assert b2.max_tokens == 8192

    resp = {
        "choices": [
            {
                "finish_reason": "tool_calls",
                "message": {
                    "role": "assistant",
                    "reasoning": "chain of thought",
                    "reasoning_content": "more chain of thought",
                    "content": None,
                    "tool_calls": [
                        {
                            "id": "call_qwen_1",
                            "type": "function",
                            "function": {"name": "get_variant_evidence", "arguments": "{}"},
                        }
                    ],
                },
            }
        ]
    }
    msg = _completion_to_message(resp)
    assert msg.content is None
    assert len(msg.tool_calls) == 1
    assert msg.tool_calls[0].id == "call_qwen_1"
    serialized = _message_to_openai(msg)
    assert "reasoning" not in serialized
    assert "reasoning_content" not in serialized
    assert serialized["tool_calls"][0]["id"] == "call_qwen_1"


def test_openai_backend_length_truncation_retries_once_2x(monkeypatch):
    from acmg_agent.llm_backends import OpenAICompatBackend

    backend = OpenAICompatBackend(
        "http://unimatrix01.kaust.edu.sa:8000/v1",
        "qwen3.8-27b",
        max_tokens=4096,
        backoff_seconds=0,
    )
    seen_max_tokens: list[int] = []
    responses = [
        {
            "choices": [
                {
                    "finish_reason": "length",
                    "message": {"role": "assistant", "reasoning": "all tokens spent reasoning", "content": None},
                }
            ]
        },
        {
            "choices": [
                {
                    "finish_reason": "tool_calls",
                    "message": {
                        "role": "assistant",
                        "reasoning": "finished thinking",
                        "content": None,
                        "tool_calls": [
                            {
                                "id": "call_1",
                                "type": "function",
                                "function": {"name": "get_variant_evidence", "arguments": "{}"},
                            }
                        ],
                    },
                }
            ]
        },
    ]

    def fake_post(body):
        seen_max_tokens.append(body["max_tokens"])
        return responses.pop(0)

    monkeypatch.setattr(backend, "_post_json", fake_post)
    reply = backend.chat([Message(role="user", content="hi")], tools=[{"type": "function"}])
    assert seen_max_tokens == [4096, 8192]
    assert len(reply.tool_calls) == 1
    assert reply.tool_calls[0].name == "get_variant_evidence"


def test_openai_backend_double_length_truncation_surfaces_failed(monkeypatch):
    from acmg_agent.llm_backends import OpenAICompatBackend

    backend = OpenAICompatBackend(
        "http://unimatrix01.kaust.edu.sa:8000/v1",
        "qwen3.8-27b",
        max_tokens=4096,
        backoff_seconds=0,
    )
    truncated = {
        "choices": [
            {
                "finish_reason": "length",
                "message": {"role": "assistant", "reasoning": "still reasoning", "content": None},
            }
        ]
    }
    monkeypatch.setattr(backend, "_post_json", lambda body: truncated)
    ev = make_ev()
    det = det_criteria(ev)
    res = run_agent(ev, det, backend)
    assert res.status == "failed"
    assert res.criteria == det


def test_openai_backend_http_retries_and_session_affinity(monkeypatch):
    import http.client
    import json as _json
    from acmg_agent.llm_backends import OpenAICompatBackend

    backend = OpenAICompatBackend(
        "http://unimatrix01.kaust.edu.sa:8000/v1",
        "qwen3.8-27b",
        sample="TS2",
        backoff_seconds=0,
    )
    assert backend.session_affinity.startswith("TS2-")
    captured_headers: list[dict] = []
    outcomes = [
        ConnectionResetError("reset"),
        (503, "{}"),
        (
            200,
            _json.dumps({
                "choices": [
                    {
                        "finish_reason": "stop",
                        "message": {"role": "assistant", "content": "ok reply"},
                    }
                ]
            }),
        ),
    ]

    class FakeResp:
        def __init__(self, status, body):
            self.status = status
            self._body = body.encode()

        def read(self):
            return self._body

    class FakeConn:
        def __init__(self, netloc, timeout=None):
            self.netloc = netloc
            self.timeout = timeout

        def request(self, method, path, body=None, headers=None):
            captured_headers.append(dict(headers or {}))
            item = outcomes.pop(0)
            if isinstance(item, Exception):
                raise item
            self._resp = FakeResp(*item)

        def getresponse(self):
            return self._resp

        def close(self):
            pass

    monkeypatch.setattr(http.client, "HTTPConnection", FakeConn)
    msg = backend.chat([Message(role="user", content="ping")], tools=[])
    assert msg.content == "ok reply"
    assert len(captured_headers) == 3
    assert all(h.get("X-Session-Affinity") == backend.session_affinity for h in captured_headers)


def test_fake_backend_narrative_finish_reason_length_retries_and_traces(tmp_path):
    import json as _json
    ev = make_ev(gene="HEXA")
    det = det_criteria(ev)
    payload = {
        "criteria": [],
        "rationale": "No adjustments were needed; PVS1, PM2, and PM3 are met.",
    }
    cut_sentence = (
        "The final criteria list for HEXA variant 1:100 A>G records PVS1 as met "
        "(very_strong) for a splice_donor_variant in a LoF-mechanism gene that is not in an NMD"
    )
    complete_sentence = (
        "The final criteria list for HEXA variant 1:100 A>G records PVS1 as met "
        "(very_strong) for a splice_donor_variant in a LoF-mechanism gene that is not in an "
        "NMD-escape region, together with PM2 and PM3 at supporting strength."
    )
    backend = FakeBackend([
        submit_tool(payload),
        Message(
            role="assistant",
            content=cut_sentence,
            finish_reason="length",
            usage={"prompt_tokens": 350, "completion_tokens": 8192, "total_tokens": 8542},
        ),
        Message(
            role="assistant",
            content=complete_sentence,
            finish_reason="stop",
            usage={"prompt_tokens": 420, "completion_tokens": 110, "total_tokens": 530},
        ),
    ])
    trace_file = tmp_path / "llm-trace.jsonl"
    backend.trace_path = str(trace_file)

    res = run_agent(ev, det, backend)
    assert res.status == "ok"
    assert res.narrative == complete_sentence
    # Initial narrative call uses >= 8192 max_tokens; retry doubles to >= 16384
    assert backend.max_tokens_seen[1] >= 8192
    assert backend.max_tokens_seen[2] >= 16384

    assert trace_file.exists()
    records = [_json.loads(line) for line in trace_file.read_text().splitlines() if line.strip()]
    assert len(records) == 3
    assert records[0]["phase"] == "review"
    assert records[1]["phase"] == "narrative"
    assert records[1]["response"]["finish_reason"] == "length"
    assert records[1]["response"]["usage"]["completion_tokens"] == 8192
    assert records[2]["phase"] == "narrative_retry"
    assert records[2]["response"]["finish_reason"] == "stop"


def test_fake_backend_narrative_double_truncation_marks_ok_truncated(tmp_path):
    from acmg_agent.report import render_markdown, serialize_candidate_assessment
    from acmg_agent.classify import classify
    import acmg_agent.pdf_report as pdf_mod

    ev = make_ev(gene="PAK2")
    det = det_criteria(ev)
    payload = {
        "criteria": [],
        "rationale": "PAK2 criteria confirmed deterministically.",
    }
    cut1 = "The PAK2 variant is annotated as a frameshift_variant. PVS1 is not met because LoF is unknown and LOEUF=0"
    cut2 = "The PAK2 variant is annotated as a frameshift_variant in the supplied criteria list and LOEUF=0"
    backend = FakeBackend([
        submit_tool(payload),
        Message(role="assistant", content=cut1, finish_reason="length"),
        Message(role="assistant", content=cut2, finish_reason="length"),
    ])
    trace_file = tmp_path / "llm-trace.jsonl"
    backend.trace_path = str(trace_file)

    res = run_agent(ev, det, backend)
    assert res.status == "ok_truncated"
    assert res.rationale == "PAK2 criteria confirmed deterministically."

    # Verify report (Markdown & PDF HTML) explicitly labels ok_truncated and does not silently print a cut sentence as ok
    clf = classify(res.criteria)
    cand_dict = serialize_candidate_assessment(
        {"gene": "PAK2", "chrom": "1", "pos": "100", "ref": "A", "alt": "G", "hgvs": "c.100+1G>A"},
        ev,
        res.criteria,
        clf,
        tier=1,
        adjusted_codes=res.adjusted,
        agent_rationale=res.rationale,
        agent_narrative=res.narrative,
        agent_status=res.status,
    )
    assert cand_dict["agent_status"] == "ok_truncated"
    report_data = {
        "schema_version": 1,
        "sample": "sample",
        "patient_id": "TS1",
        "llm": {
            "backend": "openai",
            "model": "qwen3.8-27b",
            "status": "partial",
            "reviewed": ["PAK2"],
            "unreviewed": [],
            "statuses": {"PAK2": "ok_truncated"},
        },
        "leading_candidate": cand_dict,
        "near_tie_alternatives": [],
        "candidates": [cand_dict],
        "expansions": [],
        "narrative": "Deterministic leading summary.",
        "agent_narrative": res.narrative,
        "errors": [],
    }
    md = render_markdown(report_data)
    assert "ok_truncated" in md
    assert "truncated at the token limit" in md.lower() or "truncated at token limit" in md.lower()
    assert "Patient ID**: TS1 (VCF sample: sample)" in md

    html = pdf_mod.render_pdf_html(report_data)
    assert "ok_truncated" in html
    assert "Patient ID: <code>TS1</code> (VCF sample: <code>sample</code>)" in html
    assert "Patient ID: TS1 (VCF sample: sample)" in html
