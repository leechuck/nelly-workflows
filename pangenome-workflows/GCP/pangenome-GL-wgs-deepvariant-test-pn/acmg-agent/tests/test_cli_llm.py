"""CLI -> agent wiring tests: FakeBackend injected into cmd_assess via monkeypatch.

No real model is contacted; the tests only verify that the WP4 agent layer is
invoked by the CLI, that validated adjustments surface as adjusted_by_agent
flags, and that llm.{backend,model,status} is recorded (plan Task 4.1).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import acmg_agent.cli as cli_mod
from acmg_agent.cli import main
from acmg_agent.llm_backends import Message, ToolCall

FIX = Path(__file__).parent / "fixtures" / "acmg"


class CliFakeBackend:
    """Per-gene-scripted tool-calling backend: submits once per gene, then
    answers with a content-only narrative prompt."""

    backend_name = "fake"
    model_name = "fake-model-7b"

    def __init__(self) -> None:
        self.calls: list[tuple] = []
        self.submitted_per_gene: dict[str, dict] = {}
        self._submitted: set[str] = set()

    def set_submission(self, gene: str, payload: dict) -> None:
        self.submitted_per_gene[gene] = payload

    def chat(self, messages: list[Message], tools: list[dict]) -> Message:
        self.calls.append((messages, tools))
        gene = self._gene_of(messages)
        payload = self.submitted_per_gene.get(gene)
        if payload is not None and gene not in self._submitted:
            self._submitted.add(gene)
            return Message(role="assistant", tool_calls=[ToolCall(name="submit_assessment", arguments=payload)])
        return Message(role="assistant", content=f"review narrative for {gene}")

    @staticmethod
    def _gene_of(messages: list[Message]) -> str:
        text = next((m.content or "" for m in messages if m.role == "user"), "")
        hit = re.search(r'"gene":\s*"([^"]+)"', text) or re.search(r"candidate gene ([^,\s]+)", text)
        return hit.group(1) if hit else "???"


def test_cli_runs_agent_layer_and_records_flags(tmp_path, monkeypatch):
    backend = CliFakeBackend()
    payload = {
        "criteria": [
            {"code": "PVS1", "status": "met", "strength": "strong",
             "reason": "nmd_escape per the pack"},
        ],
        "rationale": "downgrade per nmd_escape",
    }
    backend.set_submission("HEXA", payload)
    monkeypatch.setattr(cli_mod, "select_backend_from_env", lambda: backend)

    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    data = json.loads(out_json.read_text())
    # llm info recorded from the backend
    assert data["llm"]["backend"] == "fake"
    assert data["llm"]["model"] == "fake-model-7b"
    assert data["llm"]["status"] in ("ok", "partial")
    assert len(backend.calls) > 0
    hexa_items = [c for c in data["candidates"] if c["gene"] == "HEXA"]
    assert hexa_items, "HEXA must appear among candidates"
    hexa = hexa_items[0]
    pvs1 = [c for c in hexa["criteria"] if c["code"] == "PVS1"][0]
    assert pvs1["strength"] == "strong"
    assert pvs1["adjusted_by_agent"] is True
    other = [c for c in hexa["criteria"] if c["code"] != "PVS1"]
    assert all(c.get("adjusted_by_agent") is False for c in other)
    assert hexa["agent_rationale"] == "downgrade per nmd_escape"
    assert hexa["agent_status"] == "ok"
    assert hexa["adjusted_by_agent"] == ["PVS1"]
    assert data["leading_candidate"]["gene"] == "HEXA"
    assert data["leading_candidate"]["agent_status"] == "ok"
    # Defect (4): every candidate has agent_status and adjusted_by_agent, and llm has reviewed/unreviewed/statuses
    valid_statuses = {"ok", "timeout", "failed", "rejected", "skipped", "not_reviewed"}
    for c in data["candidates"]:
        assert c.get("agent_status") in valid_statuses
        assert isinstance(c.get("adjusted_by_agent"), list)
    assert data["llm"]["reviewed"] == ["HEXA"]
    assert "SEMA6B" in data["llm"]["unreviewed"]
    assert data["llm"]["statuses"]["HEXA"] == "ok"

    # Defect (1): Markdown report must include leading narrative under section 1
    # and section 5 LLM review with per-candidate status, adjusted criteria, rationale, narrative.
    md = out_md.read_text()
    assert "- **LLM Review**: backend=`fake`, model=`fake-model-7b`, status=`partial`" in md
    assert "**Leading Candidate LLM Narrative (HEXA)**: review narrative for HEXA" in md
    assert "## 5. LLM review (fake-model-7b, partial)" in md
    assert "- **Reviewed genes**: HEXA" in md
    assert "- **Status**: `ok`" in md
    assert "`PVS1` -> `met` (`strong`): nmd_escape per the pack" in md
    assert "- **Rationale**: downgrade per nmd_escape" in md
    assert "- **Narrative**: review narrative for HEXA" in md
    # First call to the backend must be for the leading candidate HEXA (Defect 3)
    first_gene = CliFakeBackend._gene_of(backend.calls[0][0])
    assert first_gene == "HEXA"


def test_cli_markdown_warns_when_leading_unreviewed(tmp_path, monkeypatch):
    backend = CliFakeBackend()
    # Only script SEMA6B; HEXA (leading) gets no submission so it is rejected/unreviewed
    payload = {
        "criteria": [
            {"code": "PP3", "status": "met", "strength": "supporting",
             "reason": "revel per the pack"},
        ],
        "rationale": "SEMA6B reviewed",
    }
    backend.set_submission("SEMA6B", payload)
    monkeypatch.setattr(cli_mod, "select_backend_from_env", lambda: backend)

    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    data = json.loads(out_json.read_text())
    assert data["leading_candidate"]["gene"] == "HEXA"
    assert data["leading_candidate"]["agent_status"] == "rejected"
    # Never mark the run 'ok' while the leading candidate is unreviewed
    assert data["llm"]["status"] == "partial"

    md = out_md.read_text()
    assert "WARNING — Leading candidate NOT reviewed by LLM" in md
    assert "Leading candidate **HEXA** was NOT reviewed by the LLM" in md
    assert "## 5. LLM review (fake-model-7b, partial)" in md
    # First candidate sent to LLM was still HEXA (leading)
    first_gene = CliFakeBackend._gene_of(backend.calls[0][0])
    assert first_gene == "HEXA"


def test_cli_no_llm_stays_disabled(tmp_path, monkeypatch):
    monkeypatch.setenv("ACMG_LLM_BACKEND", "off")
    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--no-llm",
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    data = json.loads(out_json.read_text())
    assert data["llm"]["status"] == "disabled"
    assert data["llm"]["reviewed"] == []
    assert "HEXA" in data["llm"]["unreviewed"]
    for c in data["candidates"]:
        assert c["agent_status"] == "not_reviewed"
        assert c["adjusted_by_agent"] == []


def test_cli_expansion_led_hd2_reviews_leading_expansion_and_renders_section_1_narrative(tmp_path, monkeypatch):
    backend = CliFakeBackend()
    monkeypatch.setattr(cli_mod, "select_backend_from_env", lambda: backend)

    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "HD2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "HD2/exomiser.variants.tsv"),
        "--eh-json", str(FIX / "HD2/eh_repeats.json"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    data = json.loads(out_json.read_text())
    lead = data["leading_candidate"]
    assert lead["kind"] == "expansion"
    assert lead["gene"] == "HTT"
    assert lead["genotypes"] == "18/46"
    assert lead["agent_status"] == "ok"
    assert lead["adjusted_by_agent"] == []
    assert lead["agent_narrative"] == "review narrative for HTT"
    assert data["agent_narrative"] == "review narrative for HTT"
    assert data["expansions"][0]["gene"] == "HTT"
    assert data["expansions"][0]["agent_status"] == "ok"
    assert data["expansions"][0]["agent_narrative"] == "review narrative for HTT"
    assert data["llm"]["reviewed"] == ["HTT"]
    assert data["llm"]["statuses"]["HTT"] == "ok"

    # First call to backend must be for the leading repeat expansion HTT
    first_gene = CliFakeBackend._gene_of(backend.calls[0][0])
    assert first_gene == "HTT"

    md = out_md.read_text()
    assert "**Leading Candidate LLM Narrative (HTT)**: review narrative for HTT" in md
    assert "### HTT (repeat expansion 18/46) — status: `ok`" in md
    assert "- **Adjusted criteria**: none (repeat expansion; no ACMG criteria adjustments)" in md
    assert "- **Narrative**: review narrative for HTT" in md
