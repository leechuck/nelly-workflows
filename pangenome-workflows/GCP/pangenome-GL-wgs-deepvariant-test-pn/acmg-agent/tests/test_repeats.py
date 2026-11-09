import json
from pathlib import Path
from acmg_agent.repeats import load_strchive, classify_repeat
from acmg_agent.cli import main

FIX = Path(__file__).parent / "fixtures" / "acmg"
SUBSET = FIX / "strchive-loci.subset.json"


def test_htt_46_is_pathogenic():
    loci = load_strchive(SUBSET)
    assert "HTT" in loci
    eh_var = {
        "VariantId": "HTT",
        "Genotype": "18/46",
        "GenotypeConfidenceInterval": "18-18/46-58",
        "ReferenceRegion": "chr4:3074876-3074933",
        "RepeatUnit": "CAG",
    }
    call = classify_repeat(eh_var, loci)
    assert call.gene == "HTT"
    assert call.max_allele == 46
    assert call.tier == "pathogenic"
    assert call.ranges["pathogenic"][0] == 36
    assert "Huntington" in call.disease
    assert call.moi == "AD"


def test_htt_30_is_intermediate():
    loci = load_strchive(SUBSET)
    eh_var = {
        "VariantId": "HTT",
        "Genotype": "18/30",
        "GenotypeConfidenceInterval": "18-18/30-30",
        "ReferenceRegion": "chr4:3074876-3074933",
    }
    call = classify_repeat(eh_var, loci)
    assert call.max_allele == 30
    assert call.tier == "intermediate"


def test_locus_absent_is_unknown():
    loci = load_strchive(SUBSET)
    eh_var = {
        "VariantId": "NOT_A_GENE",
        "Genotype": "10/50",
        "ReferenceRegion": "chr1:100-200",
    }
    call = classify_repeat(eh_var, loci)
    assert call.tier == "unknown_locus"

    # Also verify non-overlapping ReferenceRegion on same gene symbol (e.g. HTT_CCG) is rejected
    eh_ccg = {
        "VariantId": "HTT",
        "Genotype": "9/45",
        "ReferenceRegion": "chr4:3074939-3074966",
    }
    assert classify_repeat(eh_ccg, loci).tier == "unknown_locus"


def test_hd2_fixture_leads_with_htt_via_strchive(tmp_path):
    out_md = tmp_path / "summary.md"
    out_json = tmp_path / "acmg.json"
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "HD2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "HD2/exomiser.variants.tsv"),
        "--eh-json", str(FIX / "HD2/eh_repeats.json"),
        "--strchive-json", str(SUBSET),
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
    assert data["leading_candidate"]["allele2"] == 46
    assert data["leading_candidate"]["status"] == "pathogenic_range"
    assert data.get("strchive_used") is True
