import json
from pathlib import Path
import pysam
from evaluate_erepo import (
    build_erepo_vcf,
    compare_assessment_to_erepo,
    main as erepo_main,
)

FIX = Path(__file__).parent / "fixtures" / "acmg"
EREPO_SUBSET = FIX / "erepo.subset.tsv"


def test_erepo_vcf_builder_parses_nc_hgvs(tmp_path):
    out_vcf = tmp_path / "erepo.vcf.gz"
    rc = erepo_main([
        "--erepo", str(EREPO_SUBSET),
        "--out-vcf", str(out_vcf),
    ])
    assert rc == 0
    assert out_vcf.exists()
    with pysam.VariantFile(str(out_vcf)) as vf:
        recs = list(vf)
    # Out of 200 rows in erepo.subset.tsv, >150 have GRCh38 NC_ SNV expressions
    assert len(recs) >= 150
    first_positions = {(r.chrom, r.pos, r.ref, r.alts[0]) for r in recs}
    # Row 2 in erepo.subset.tsv: NC_000012.12:g.102917130T>C
    assert ("chr12", 102917130, "T", "C") in first_positions


def test_concordance_table_counts(tmp_path):
    tsv_path = tmp_path / "mini_erepo.tsv"
    tsv_path.write_text(
        "#Variation\tClinVar Variation Id\tAllele Registry Id\tHGVS Expressions\tHGNC Gene Symbol\tDisease\tMondo Id\tMode of Inheritance\tAssertion\tApplied Evidence Codes (Met)\tApplied Evidence Codes (Not Met)\n"
        "var1\t1001\tCA1\tNC_000012.12:g.102917130T>C\tPAH\tPKU\tMONDO:1\tAutosomal recessive inheritance\tPathogenic\tPVS1, PM2_Supporting, PP3_Moderate\tBS1, BA1\n"
        "var2\t1002\tCA2\tNC_000012.12:g.102840437A>G\tPAH\tPKU\tMONDO:1\tAutosomal recessive inheritance\tBenign\tBA1, BP4\tPM2\n"
    )
    assessment_path = tmp_path / "acmg-assessment.json"
    assessment_data = {
        "schema_version": 1,
        "candidates": [
            {
                "gene": "PAH",
                "variant": "chr12:102917130 T>C",
                "classification": "Pathogenic",
                "points": 11,
                "criteria": [
                    {"code": "PVS1", "status": "met", "strength": "very_strong", "rationale": "r"},
                    {"code": "PM2", "status": "met", "strength": "supporting", "rationale": "r"},
                    {"code": "PP3", "status": "not_met", "strength": None, "rationale": "r"},
                    {"code": "BS1", "status": "not_met", "strength": None, "rationale": "r"},
                    {"code": "BA1", "status": "not_met", "strength": None, "rationale": "r"},
                ],
            },
            {
                "gene": "PAH",
                "variant": "chr12:102840437 A>G",
                "classification": "Benign",
                "points": -9,
                "criteria": [
                    {"code": "BA1", "status": "met", "strength": "stand_alone", "rationale": "r"},
                    {"code": "BP4", "status": "met", "strength": "supporting", "rationale": "r"},
                    {"code": "PM2", "status": "not_met", "strength": None, "rationale": "r"},
                ],
            },
        ],
    }
    assessment_path.write_text(json.dumps(assessment_data))
    out_md = tmp_path / "erepo-concordance.md"
    metrics = compare_assessment_to_erepo(assessment_path, tsv_path, out_md)
    assert out_md.exists()
    # PVS1: 1 TP, 0 FP, 0 FN
    assert metrics["per_criterion"]["PVS1"]["tp"] == 1
    assert metrics["per_criterion"]["PVS1"]["fn"] == 0
    # PM2 (normalized from PM2_Supporting): 1 TP, 0 FP, 0 FN, 1 TN
    assert metrics["per_criterion"]["PM2"]["tp"] == 1
    assert metrics["per_criterion"]["PM2"]["tn"] == 1
    # PP3 (normalized from PP3_Moderate): engine said not_met -> 1 FN
    assert metrics["per_criterion"]["PP3"]["fn"] == 1
    # BA1: 1 TP (var2), 1 TN (var1)
    assert metrics["per_criterion"]["BA1"]["tp"] == 1
    assert metrics["per_criterion"]["BA1"]["tn"] == 1
    # Class concordance: 2/2 exact matches
    assert metrics["class_matrix"][("Pathogenic", "Pathogenic")] == 1
    assert metrics["class_matrix"][("Benign", "Benign")] == 1
