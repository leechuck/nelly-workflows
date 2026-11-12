"""R12 Task 1: the LLM reviewer must not recommend clinical actions."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import acmg_agent.cli as cli_mod
from acmg_agent.cli import main
from acmg_agent.llm_backends import Message
from acmg_agent.policy import REDACTION_NOTE, find_recommendations, strip_recommendations
from tests.test_cli_llm import CliFakeBackend
from tests.test_pdf_report import _extract_pdf_text

FIX = Path(__file__).parent / "fixtures" / "acmg"

# Verbatim from the real Huntington run (Qwen3.8-27B, v7, 2026-11-11).
HD_REAL = ("The result supports a finding consistent with Huntington disease and is amenable to "
           "genetic counseling regarding predictive testing for at-risk family members.")
BAD = "Genetic counselling and predictive testing of relatives should be offered."


def test_real_hd_sentence_is_caught():
    assert find_recommendations(HD_REAL) == [HD_REAL]


@pytest.mark.parametrize("text", [
    "The allele is in the STRchive full-penetrance range.",
    "PVS1 applies because the splice donor variant is in a gene where loss of function is the disease mechanism.",
    "PM3_Supporting is met for a homozygous variant in an autosomal recessive gene.",
    "No treatment-related evidence is in the pack.",
])
def test_neutral_sentences_pass(text):
    assert find_recommendations(text) == []


@pytest.mark.parametrize("text", [
    BAD,
    "We recommend confirmatory testing.",
    "Cascade testing of relatives is advised.",
    "Referral to a neurologist should be considered.",
    "Consider carrier testing for the partner.",
    "Clinical management should include surveillance.",
])
def test_recommendations_caught(text):
    assert find_recommendations(text)


def test_strip_keeps_neutral_sentences():
    clean, removed = strip_recommendations(f"The HTT alleles are 18/46. {BAD}")
    assert removed and clean == "The HTT alleles are 18/46."


class PolicyBackend(CliFakeBackend):
    """Narrative replies append a recommendation; the policy-retry reply is clean or not."""

    def __init__(self, clean_on_retry: bool):
        super().__init__()
        self.clean_on_retry = clean_on_retry
        self.violation_turns = 0

    def chat(self, messages: list[Message], tools: list[dict], **_kw) -> Message:
        retry = any("POLICY VIOLATION" in (m.content or "") for m in messages)
        msg = super().chat(messages, tools)
        if msg.tool_calls:
            return msg
        if retry:
            self.violation_turns += 1
            gene = "HTT" if "HTT" in json.dumps([m.content for m in messages]) else "HEXA"
            text = f"Neutral narrative for {gene}." if self.clean_on_retry else f"Rewritten narrative for {gene}. {BAD}"
            return Message(role="assistant", content=text, finish_reason="stop")
        return Message(role="assistant", content=f"{msg.content} {BAD}", finish_reason="stop")


def _run(tmp_path, monkeypatch, backend, case):
    backend.set_submission("HEXA", {"criteria": [], "rationale": "No change: the pack supports the engine criteria."})
    monkeypatch.setattr(cli_mod, "select_backend_from_env", lambda: backend)
    args = ["assess",
            "--genes-tsv", str(FIX / f"{case}/exomiser.genes.tsv"),
            "--variants-tsv", str(FIX / f"{case}/exomiser.variants.tsv"),
            "--vep-vcf", str(FIX / "vep.vcf"),
            "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
            "--gene-table", str(FIX / "gene_table.tsv"),
            "--output-md", str(tmp_path / "r.md"),
            "--output-json", str(tmp_path / "a.json"),
            "--output-pdf", str(tmp_path / "r.pdf")]
    if case == "HD2":
        args[1:1] = ["--eh-json", str(FIX / "HD2/eh_repeats.json")]
    assert main(args) == 0
    data = json.loads((tmp_path / "a.json").read_text())
    return data, (tmp_path / "r.md").read_text(), _extract_pdf_text(tmp_path / "r.pdf")[1]


def _all_llm_text(data):
    items = list(data["candidates"]) + list(data.get("expansions") or [])
    return " ".join(str(i.get(k) or "") for i in items for k in ("agent_narrative", "agent_rationale"))


@pytest.mark.parametrize("case,gene", [("TS2", "HEXA"), ("HD2", "HTT")])
def test_retry_gives_clean_text(tmp_path, monkeypatch, case, gene):
    data, md, pdf = _run(tmp_path, monkeypatch, PolicyBackend(clean_on_retry=True), case)
    lead = data["leading_candidate"]
    assert lead["gene"] == gene
    assert lead["agent_narrative"] == f"Neutral narrative for {gene}."
    assert not lead.get("narrative_redacted")
    assert REDACTION_NOTE not in md


@pytest.mark.parametrize("case,gene", [("TS2", "HEXA"), ("HD2", "HTT")])
def test_persistent_violation_is_redacted_and_noted(tmp_path, monkeypatch, case, gene):
    data, md, pdf = _run(tmp_path, monkeypatch, PolicyBackend(clean_on_retry=False), case)
    lead = data["leading_candidate"]
    assert lead["gene"] == gene
    assert lead.get("narrative_redacted") is True
    assert lead["agent_narrative"] and gene in lead["agent_narrative"]
    assert find_recommendations(_all_llm_text(data)) == []
    assert "should be offered" not in md and "should be offered" not in pdf
    assert REDACTION_NOTE in md
    assert "clinical recommendations were removed" in pdf
