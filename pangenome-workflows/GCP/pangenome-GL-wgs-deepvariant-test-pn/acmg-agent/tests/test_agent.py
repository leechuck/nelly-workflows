"""Tests for the bounded LLM agent layer (FakeBackend only -- no real model)."""
from __future__ import annotations

from typing import Any

from acmg_agent.agent import (
    run_agent,
    run_all,
    run_llm_layer,
    select_review_set,
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

    Script items: ToolCall (one call), str (plain assistant content),
    Exception instance (raised). When the script is exhausted, a canned
    content-only reply is returned (narrative / finalisation nudge).
    """

    def __init__(self, script):
        self.script = list(script)
        self.calls: list[tuple[list[Message], list[dict]]] = []

    def chat(self, messages: list[Message], tools: list[dict]) -> Message:
        self.calls.append((messages, tools))
        if not self.script:
            return Message(role="assistant", content=f"narrative for candidate {len(self.calls)}")
        item = self.script.pop(0)
        if isinstance(item, Exception):
            raise item
        if isinstance(item, str):
            return Message(role="assistant", content=item)
        return Message(role="assistant", tool_calls=[item])


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
