"""Tests for clinical PDF report generation (WeasyPrint + pypdf, Task R8 / design §11)."""
from __future__ import annotations

from pathlib import Path
import pypdf

import acmg_agent.cli as cli_mod
from acmg_agent.cli import main
import acmg_agent.pdf_report as pdf_mod
from tests.test_cli_llm import CliFakeBackend

FIX = Path(__file__).parent / "fixtures" / "acmg"


def _extract_pdf_text(pdf_path: Path) -> tuple[int, str]:
    reader = pypdf.PdfReader(str(pdf_path))
    pages = len(reader.pages)
    text = "\n".join(page.extract_text() or "" for page in reader.pages)
    return pages, text


def test_pdf_report_ts2_deterministic(tmp_path):
    out_md = tmp_path / "acmg-report.md"
    out_json = tmp_path / "acmg-assessment.json"
    out_pdf = tmp_path / "acmg-report.pdf"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--sample", "TS2",
        "--sex", "FEMALE",
        "--hpo-terms", "HP:0001263,HP:0001250",
        "--phenotype", "Global developmental delay; Seizures",
        "--no-llm",
        "--output-md", str(out_md),
        "--output-json", str(out_json),
        "--output-pdf", str(out_pdf),
    ])
    assert rc == 0
    assert out_pdf.exists() and out_pdf.stat().st_size > 1000
    n_pages, text = _extract_pdf_text(out_pdf)
    assert n_pages >= 2
    for expected in ("PVS1", "Met", "HEXA", "Appendix"):
        assert expected in text, f"Missing {expected!r} in TS2 PDF text"


def test_pdf_report_hd2_deterministic(tmp_path):
    out_md = tmp_path / "acmg-report.md"
    out_json = tmp_path / "acmg-assessment.json"
    out_pdf = tmp_path / "acmg-report.pdf"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "HD2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "HD2/exomiser.variants.tsv"),
        "--eh-json", str(FIX / "HD2/eh_repeats.json"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--sample", "HD2",
        "--sex", "MALE",
        "--hpo-terms", "HP:0002072,HP:0000726",
        "--phenotype", "Chorea; Dementia",
        "--no-llm",
        "--output-md", str(out_md),
        "--output-json", str(out_json),
        "--output-pdf", str(out_pdf),
    ])
    assert rc == 0
    assert out_pdf.exists() and out_pdf.stat().st_size > 1000
    n_pages, text = _extract_pdf_text(out_pdf)
    assert n_pages >= 2
    for expected in ("PVS1", "Met", "HTT", "Appendix"):
        assert expected in text, f"Missing {expected!r} in HD2 PDF text"


def test_pdf_report_ts2_with_fake_backend_rationale(tmp_path, monkeypatch):
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

    out_md = tmp_path / "acmg-report.md"
    out_json = tmp_path / "acmg-assessment.json"
    out_pdf = tmp_path / "acmg-report.pdf"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--output-md", str(out_md),
        "--output-json", str(out_json),
        "--output-pdf", str(out_pdf),
    ])
    assert rc == 0
    assert out_pdf.exists()
    n_pages, text = _extract_pdf_text(out_pdf)
    assert n_pages >= 2
    for expected in ("PVS1", "Met", "HEXA", "Appendix", "downgrade per nmd_escape", "review narrative for HEXA"):
        assert expected in text, f"Missing {expected!r} in TS2 LLM PDF text"


def test_pdf_failure_writes_error_note_and_exits_zero(tmp_path, monkeypatch):
    def boom(path, data):
        raise RuntimeError("simulated WeasyPrint failure")

    monkeypatch.setattr(pdf_mod, "write_pdf", boom)

    out_md = tmp_path / "acmg-report.md"
    out_json = tmp_path / "acmg-assessment.json"
    out_pdf = tmp_path / "acmg-report.pdf"
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
        "--output-pdf", str(out_pdf),
    ])
    assert rc == 0
    assert out_json.exists()
    assert out_md.exists()
    err_file = Path(f"{out_pdf}.error.txt")
    assert err_file.exists()
    assert "simulated WeasyPrint failure" in err_file.read_text()
