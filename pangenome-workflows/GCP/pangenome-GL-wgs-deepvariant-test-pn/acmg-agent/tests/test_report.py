from pathlib import Path
import json
from acmg_agent.cli import main

FIX = Path(__file__).parent / "fixtures" / "acmg"


def test_cli_no_llm_ts2_hexa_leading(tmp_path):
    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--eh-json", str(FIX / "TS2/eh_repeats.json"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--no-llm",
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    data = json.loads(out_json.read_text())
    assert data["schema_version"] == 1
    assert data["leading_candidate"]["gene"] == "HEXA"
    assert data["leading_candidate"]["classification"] == "Pathogenic"
    assert data["leading_candidate"]["points"] == 12
    md = out_md.read_text()
    assert "HEXA" in md
    assert "PVS1" in md
    assert "PM3" in md


def test_cli_hd2_htt_expansion_leading(tmp_path):
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
        "--no-llm",
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    data = json.loads(out_json.read_text())
    assert data["leading_candidate"]["kind"] == "expansion"
    assert data["leading_candidate"]["gene"] == "HTT"


def test_outputs_always_written_on_bad_vep(tmp_path):
    bad_vep = tmp_path / "broken.vcf"
    bad_vep.write_text("not a vcf")
    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--vep-vcf", str(bad_vep),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--no-llm",
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    assert out_md.exists() and out_json.exists()


def test_extract_writes_candidates_vcf(tmp_path):
    out_vcf = tmp_path / "cands.vcf"
    rc = main([
        "extract",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--fai", str(FIX / "hg38.fa.fai"),
        "--output-vcf", str(out_vcf),
    ])
    assert rc == 0
    text = out_vcf.read_text()
    assert "72346234" in text


def test_from_vcf_mode_no_exomiser(tmp_path):
    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--from-vcf",
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--no-llm",
        "--output-md", str(out_md),
        "--output-json", str(out_json),
    ])
    assert rc == 0
    data = json.loads(out_json.read_text())
    assert len(data["candidates"]) >= 1
