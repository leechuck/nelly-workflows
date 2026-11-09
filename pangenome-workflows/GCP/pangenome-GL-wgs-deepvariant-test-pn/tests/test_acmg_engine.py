"""End-to-end deterministic ACMG engine tests on TS, TS2, HD, and HD2 fixtures."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "acmg-agent"))

from acmg_agent.cli import main  # noqa: E402

FIX_CASE = HERE / "fixtures"
FIX_ACMG = HERE.parent / "acmg-agent" / "tests" / "fixtures" / "acmg"


def _run_case(case_name: str, tmp_path: Path) -> tuple[dict, str]:
    out_md = tmp_path / f"{case_name}.summary.md"
    out_json = tmp_path / f"{case_name}.acmg.json"
    case_dir = FIX_CASE / case_name
    rc = main([
        "assess",
        "--genes-tsv", str(case_dir / "exomiser.genes.tsv"),
        "--variants-tsv", str(case_dir / "exomiser.variants.tsv"),
        "--eh-json", str(case_dir / "eh_repeats.json"),
        "--vep-vcf", str(FIX_ACMG / "vep.vcf"),
        "--clinvar-vcf", str(FIX_ACMG / "clinvar.vcf.gz"),
        "--gene-table", str(FIX_ACMG / "gene_table.tsv"),
        "--slivar-vcf", str(FIX_ACMG / "slivar.vcf"),
        "--no-llm",
        "--sample", case_name,
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    assert out_md.exists() and out_json.exists()
    return json.loads(out_json.read_text()), out_md.read_text()


def test_ts_end_to_end(tmp_path):
    data, md = _run_case("TS", tmp_path)
    assert data["schema_version"] == 1
    lead = data["leading_candidate"]
    assert lead["kind"] == "snv"
    assert lead["gene"] == "HEXA"
    assert lead["classification"] == "Pathogenic"
    assert lead["points"] == 11
    assert "HEXA" in md


def test_ts2_end_to_end(tmp_path):
    data, md = _run_case("TS2", tmp_path)
    assert data["schema_version"] == 1
    lead = data["leading_candidate"]
    assert lead["kind"] == "snv"
    assert lead["gene"] == "HEXA"
    assert lead["classification"] == "Pathogenic"
    assert lead["points"] == 10
    assert lead["rule_based_class"] == "Pathogenic"
    assert lead["conflict"] is False
    hexa_entry = next(c for c in data["candidates"] if c["gene"] == "HEXA")
    assert len(hexa_entry["criteria"]) == 28
    by_code = {c["code"]: c for c in hexa_entry["criteria"]}
    assert by_code["PVS1"]["status"] == "met" and by_code["PVS1"]["strength"] == "very_strong"
    assert by_code["PM2"]["status"] == "met" and by_code["PM2"]["strength"] == "supporting"
    assert by_code["PM3"]["status"] == "met" and by_code["PM3"]["strength"] == "supporting"
    assert by_code["PP4"]["status"] == "not_met"
    assert by_code["PP5"]["status"] == "not_applicable"
    assert by_code["BP6"]["status"] == "not_applicable"
    assert "PVS1" in md and "PM3" in md


def test_hd_end_to_end(tmp_path):
    data, md = _run_case("HD", tmp_path)
    assert data["schema_version"] == 1
    lead = data["leading_candidate"]
    assert lead["kind"] == "expansion"
    assert lead["gene"] == "HTT"
    assert lead["allele2"] == 46
    assert lead["status"] == "pathogenic_range"
    assert "HTT" in md


def test_hd2_end_to_end(tmp_path):
    data, md = _run_case("HD2", tmp_path)
    assert data["schema_version"] == 1
    lead = data["leading_candidate"]
    assert lead["kind"] == "expansion"
    assert lead["gene"] == "HTT"
    assert lead["allele2"] == 46
    assert lead["status"] == "pathogenic_range"
    assert "HTT" in md
