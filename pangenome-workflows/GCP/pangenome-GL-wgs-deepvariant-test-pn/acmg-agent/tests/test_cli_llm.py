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
        hit = re.search(r'"gene":\s*"([^"]+)"', text)
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
