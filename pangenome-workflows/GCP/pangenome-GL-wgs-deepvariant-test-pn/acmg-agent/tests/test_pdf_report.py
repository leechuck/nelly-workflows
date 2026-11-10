"""Tests for clinical PDF report generation (WeasyPrint + pypdf, Task R8 / R9 / design §11)."""
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


def test_pdf_report_ts2_deterministic(tmp_path, monkeypatch):
    monkeypatch.setenv("ACMG_PIPELINE_VERSION", "nelly-workflows modernise-2026-11")
    monkeypatch.setenv("ACMG_IMAGE_TAG", "nelly-acmg-agent:v5")
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

    # R9 (3): Title is unnumbered; sections start at 1 = Clinical Conclusion
    assert "1. Clinical Variant Interpretation Report" not in text
    assert "Clinical Variant Interpretation Report" in text
    assert "1. Clinical Conclusion" in text
    assert "2. ACMG/AMP 28-Criteria Evaluation" in text
    assert "3. Ranked SNV / Indel Candidates" in text
    assert "4. Short Tandem Repeat Expansions" in text

    # R9 (1): No raw Markdown in Summary
    assert "**HEXA**" not in text
    assert "`chr15" not in text

    # R9 (2): Absent/zero-AF variant says 'not observed in gnomAD v4 (exomes/genomes)', never 0.00e+00
    assert "not observed in gnomAD v4 (exomes/genomes)" in text
    assert "0.00e+00" not in text

    # R9 (4): TS2 has second P/LP candidate (SEMA6B Likely pathogenic) -> compact table in section 3 + Appendix C
    assert "Appendix C" in text
    assert "SEMA6B" in text

    # R9 (6): Pipeline version and image tag from env
    assert "nelly-workflows modernise-2026-11" in text
    assert "nelly-acmg-agent:v5" in text

    for expected in ("PVS1", "Met", "HEXA", "Appendix A", "Appendix B"):
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
    assert "1. Clinical Variant Interpretation Report" not in text
    assert "1. Clinical Conclusion" in text
    assert "**HTT**" not in text
    assert "not observed in gnomAD v4 (exomes/genomes)" in text
    assert "0.00e+00" not in text
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
    # R9 (7): Appendix A summary status table headers + rationale/narrative
    for expected in ("PVS1", "Met", "HEXA", "Appendix", "Wall Time", "downgrade per nmd_escape", "review narrative for HEXA"):
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
