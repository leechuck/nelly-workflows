from pathlib import Path
import pytest
from acmg_agent.model import Evidence
from acmg_agent.criteria import (
    evaluate_all, eval_pvs1, eval_ps1, eval_pm1, eval_pm2, eval_pm3,
    eval_pm4, eval_pm5, eval_pp2, eval_pp3, eval_pp4, eval_pp5,
    eval_ba1, eval_bs1, eval_bs2, eval_bp1, eval_bp2, eval_bp3,
    eval_bp4, eval_bp5, eval_bp6, eval_bp7,
)


def make_ev(**overrides) -> Evidence:
    base = dict(
        gene="TEST", moi="AR", chrom="1", pos=100, ref="A", alt="G",
        zygosity="het", transcript="NM_000001.1", mane_select=True,
        consequence="missense_variant", hgvsc="c.10A>G", hgvsp="p.Lys4Glu",
        exon="2/10", intron=None, cds_position=10, protein_position=4,
        nmd_escape=False, in_repeat_region=False,
        max_af=0.0, popmax_af=0.0, nhomalt=0, bs1_threshold=0.005,
        revel=0.5, spliceai=0.05,
        clinvar_id=None, clinvar_sig=None, clinvar_stars=0,
        same_aa_ClinVar_PLP=0, diff_aa_ClinVar_PLP=0,
        lof_mechanism=True, loeuf=0.2, pli=0.95,
        clinvar_missense_plp=5, clinvar_missense_blb=10, in_hotspot=False,
        other_var_in_gene=[], human_pheno_score=0.3, pheno_evidence="", disease="",
    )
    base.update(overrides)
    return Evidence(**base)


def test_pvs1_met_on_canonical_splice_donor_in_lof_gene():
    c = eval_pvs1(make_ev(consequence="splice_donor_variant", lof_mechanism=True, nmd_escape=False))
    assert (c.status, c.strength) == ("met", "very_strong")


def test_pvs1_downgraded_on_nmd_escape():
    c = eval_pvs1(make_ev(consequence="stop_gained", lof_mechanism=True, nmd_escape=True))
    assert (c.status, c.strength) == ("met", "strong")


def test_pvs1_last_exon_lt_10pct_is_moderate():
    c = eval_pvs1(make_ev(
        consequence="stop_gained",
        lof_mechanism=True,
        exon="12/12",
        protein_position=480,
        protein_length=500,
        nmd_escape=True,
    ))
    assert (c.status, c.strength) == ("met", "moderate")
    c2 = eval_pvs1(make_ev(
        consequence="stop_gained",
        lof_mechanism=True,
        exon="4/4",
        protein_position=100,
        protein_length=200,
        nmd_escape=True,
    ))
    assert (c2.status, c2.strength) == ("met", "strong")


def test_pvs1_not_met_on_missense():
    assert eval_pvs1(make_ev(consequence="missense_variant")).status == "not_met"


def test_pvs1_not_evaluable_when_mechanism_and_loeuf_unknown():
    assert eval_pvs1(make_ev(consequence="stop_gained", lof_mechanism=None, loeuf=None)).status == "not_evaluable"


def test_pvs1_loeuf_fallback_when_mechanism_unknown():
    c = eval_pvs1(make_ev(consequence="stop_gained", lof_mechanism=None, loeuf=0.20, nmd_escape=False))
    assert (c.status, c.strength) == ("met", "very_strong")
    assert "LOEUF" in c.rationale
    c_high = eval_pvs1(make_ev(consequence="stop_gained", lof_mechanism=None, loeuf=0.60, nmd_escape=False))
    assert c_high.status == "not_met"


def test_ps1_met_when_same_aa_2star():
    assert eval_ps1(make_ev(same_aa_ClinVar_PLP=2)).strength == "strong"
    assert eval_ps1(make_ev(same_aa_ClinVar_PLP=1)).strength == "moderate"
    assert eval_ps1(make_ev(same_aa_ClinVar_PLP=0)).status == "not_met"


def test_pm5_met_when_diff_aa_2star():
    assert eval_pm5(make_ev(diff_aa_ClinVar_PLP=2)).strength == "moderate"
    assert eval_pm5(make_ev(diff_aa_ClinVar_PLP=1)).strength == "supporting"


def test_pm1_hotspot():
    assert eval_pm1(make_ev(in_hotspot=True)).status == "met"
    assert eval_pm1(make_ev(in_hotspot=False)).status == "not_met"
    assert eval_pm1(make_ev(in_hotspot=None)).status == "not_evaluable"


def test_pm2_supporting_under_threshold():
    c = eval_pm2(make_ev(moi="AR", max_af=5e-5))
    assert (c.status, c.strength) == ("met", "supporting")


def test_pm2_not_met_common_in_ad():
    assert eval_pm2(make_ev(moi="AD", max_af=5e-5)).status == "not_met"


def test_pm3_hom_in_ar_is_supporting():
    c = eval_pm3(make_ev(moi="AR", zygosity="hom", max_af=1e-5))
    assert (c.status, c.strength) == ("met", "supporting")
    assert "0.5 pts" in c.rationale


def test_pm3_point_table_variants():
    c_trans_p = eval_pm3(make_ev(moi="AR", zygosity="het", other_var_in_gene=[{"clinvar_sig": "Pathogenic", "phased_in_trans": True}]))
    assert (c_trans_p.status, c_trans_p.strength) == ("met", "moderate")
    c_unph_p = eval_pm3(make_ev(moi="AR", zygosity="het", other_var_in_gene=[{"clinvar_sig": "Pathogenic", "gt": "0/1"}]))
    assert (c_unph_p.status, c_unph_p.strength) == ("met", "supporting")
    c_trans_lp = eval_pm3(make_ev(moi="AR", zygosity="het", other_var_in_gene=[{"clinvar_sig": "Likely_pathogenic", "phased_in_trans": True}]))
    assert (c_trans_lp.status, c_trans_lp.strength) == ("met", "supporting")
    c_unph_lp = eval_pm3(make_ev(moi="AR", zygosity="het", other_var_in_gene=[{"clinvar_sig": "Likely_pathogenic", "gt": "0/1"}]))
    assert c_unph_lp.status == "not_met"


def test_pm3_not_applicable_in_ad():
    assert eval_pm3(make_ev(moi="AD", zygosity="hom")).status == "not_met"


def test_pm4_inframe_outside_repeat():
    assert eval_pm4(make_ev(consequence="inframe_deletion", in_repeat_region=False)).status == "met"
    assert eval_pm4(make_ev(consequence="inframe_deletion", in_repeat_region=True)).status == "not_met"


def test_bp3_inframe_in_repeat():
    assert eval_bp3(make_ev(consequence="inframe_deletion", in_repeat_region=True)).status == "met"


def test_pp2_met_when_missense_enriched():
    assert eval_pp2(make_ev(clinvar_missense_plp=30, clinvar_missense_blb=3)).status == "met"
    assert eval_pp2(make_ev(clinvar_missense_plp=5, clinvar_missense_blb=10)).status == "not_met"


def test_pp3_graded_revel_and_alphamissense_and_spliceai():
    assert eval_pp3(make_ev(revel=0.95, spliceai=0.05)).strength == "strong"
    assert eval_pp3(make_ev(revel=0.80, spliceai=0.05)).strength == "moderate"
    assert eval_pp3(make_ev(revel=0.68, spliceai=0.05)).strength == "supporting"
    assert eval_pp3(make_ev(revel=0.50, spliceai=0.05)).status == "not_met"
    # SpliceAI >= 0.2 -> supporting
    c_sp = eval_pp3(make_ev(revel=None, spliceai=0.25))
    assert (c_sp.status, c_sp.strength) == ("met", "supporting")
    # AlphaMissense fallback when REVEL is None
    c_am = eval_pp3(make_ev(revel=None, alphamissense=0.65, spliceai=0.01))
    assert (c_am.status, c_am.strength) == ("met", "supporting")
    # Stronger of REVEL and SpliceAI wins without stacking
    c_both = eval_pp3(make_ev(revel=0.95, spliceai=0.80))
    assert (c_both.status, c_both.strength) == ("met", "strong")


def test_bp4_graded_revel_and_alphamissense_and_spliceai():
    assert eval_bp4(make_ev(revel=0.002, spliceai=0.02)).strength == "strong"
    assert eval_bp4(make_ev(revel=0.010, spliceai=0.02)).strength == "moderate"
    assert eval_bp4(make_ev(revel=0.150, spliceai=0.02)).strength == "supporting"
    assert eval_bp4(make_ev(revel=0.250, spliceai=0.02)).status == "not_met"
    # AlphaMissense fallback when REVEL is None
    c_am = eval_bp4(make_ev(revel=None, alphamissense=0.20, spliceai=0.02))
    assert (c_am.status, c_am.strength) == ("met", "supporting")
    # High SpliceAI blocks BP4 even when REVEL is low
    assert eval_bp4(make_ev(revel=0.002, spliceai=0.50)).status == "not_met"


def test_pp4_met_with_specific_phenotype():
    assert eval_pp4(make_ev(human_pheno_score=0.65, disease="Tay-Sachs (OMIM:272800)")).status == "met"
    assert eval_pp4(make_ev(human_pheno_score=0.55, disease="Tay-Sachs (OMIM:272800)", pheno_evidence="OMIM:272800")).status == "not_met"
    assert eval_pp4(make_ev(human_pheno_score=0.3, disease="Tay-Sachs")).status == "not_met"


def test_pp5_and_bp6_always_not_applicable():
    assert eval_pp5(make_ev(clinvar_sig="Pathogenic", clinvar_stars=4)).status == "not_applicable"
    assert eval_bp6(make_ev(clinvar_sig="Benign", clinvar_stars=4)).status == "not_applicable"


def test_ba1_met():
    assert eval_ba1(make_ev(max_af=0.06)).strength == "stand_alone"
    assert eval_ba1(make_ev(max_af=0.01)).status == "not_met"


def test_bs1_uses_gene_threshold():
    assert eval_bs1(make_ev(max_af=0.008, bs1_threshold=0.005)).status == "met"
    assert eval_bs1(make_ev(max_af=0.002, bs1_threshold=0.005)).status == "not_met"


def test_bs2_ad_one_hom():
    assert eval_bs2(make_ev(moi="AD", nhomalt=1)).status == "met"
    assert eval_bs2(make_ev(moi="AR", nhomalt=1)).status == "not_met"
    assert eval_bs2(make_ev(moi="AR", nhomalt=3)).status == "met"


def test_bp1_missense_in_truncating_only_gene():
    assert eval_bp1(make_ev(lof_mechanism=True, clinvar_missense_plp=1, clinvar_lof_plp=50)).status == "met"
    assert eval_bp1(make_ev(lof_mechanism=True, clinvar_missense_plp=30, clinvar_lof_plp=50)).status == "not_met"


def test_bp2_second_plp_in_dominant_gene():
    other = [{"clinvar_sig": "Pathogenic", "gt": "0/1"}]
    assert eval_bp2(make_ev(moi="AD", other_var_in_gene=other)).status == "met"
    assert eval_bp2(make_ev(moi="AR", other_var_in_gene=other)).status == "not_met"


def test_bp5_default_not_met():
    assert eval_bp5(make_ev()).status == "not_met"


def test_bp7_synonymous_low_spliceai():
    assert eval_bp7(make_ev(consequence="synonymous_variant", spliceai=0.01)).status == "met"
    assert eval_bp7(make_ev(consequence="missense_variant", spliceai=0.01)).status == "not_met"


def test_always_not_evaluable_codes():
    ev = make_ev()
    by_code = {c.code: c for c in evaluate_all(ev)}
    for code in ("PS2", "PS3", "PS4", "PM6", "PP1", "BS3", "BS4"):
        assert by_code[code].status == "not_evaluable"


def test_hexa_criteria():
    from acmg_agent.inputs import VepIndex, ClinvarIndex, load_gene_table
    from acmg_agent.evidence import build_evidence
    FIX = Path(__file__).parent / "fixtures" / "acmg"
    cand = {
        "gene": "HEXA", "moi": "AR", "chrom": "15", "pos": "72346234",
        "ref": "C", "alt": "G", "gt": "1/1", "effect": "splice_donor_variant",
        "hgvs": "HEXA:NM_000520.6:c.1421+1G>C:p.?", "rsid": "rs147324677",
        "clinvar": "PATHOGENIC", "max_freq": "6.3013897E-4",
        "path_score": "0.946", "pheno_score": "0.9815", "combined": "0.9994",
        "disease": "Tay-Sachs disease (OMIM:272800)", "pheno_evidence": "",
    }
    ev = build_evidence(cand, [], VepIndex.from_vcf(FIX / "vep.vcf"),
                        ClinvarIndex(FIX / "clinvar.vcf.gz"),
                        load_gene_table(FIX / "gene_table.tsv"))
    by_code = {c.code: c for c in evaluate_all(ev)}
    assert by_code["PVS1"].status == "met" and by_code["PVS1"].strength == "very_strong"
    assert by_code["PM2"].status == "met" and by_code["PM2"].strength == "supporting"
    assert by_code["PM3"].status == "met" and by_code["PM3"].strength == "supporting"
    assert by_code["PP4"].status == "met" and by_code["PP4"].strength == "supporting"
    assert by_code["PP5"].status == "not_applicable"
    assert by_code["BP6"].status == "not_applicable"
    for b in ("BA1", "BS1", "BS2", "BP4", "BP7"):
        assert by_code[b].status == "not_met"


def test_returns_28_in_order():
    codes = [c.code for c in evaluate_all(make_ev())]
    assert codes == [
        "PVS1","PS1","PS2","PS3","PS4",
        "PM1","PM2","PM3","PM4","PM5","PM6",
        "PP1","PP2","PP3","PP4","PP5",
        "BA1","BS1","BS2","BS3","BS4",
        "BP1","BP2","BP3","BP4","BP5","BP6","BP7",
    ]


def test_all_none_evidence_never_raises():
    ev = make_ev(
        transcript=None, mane_select=False, hgvsc=None, hgvsp=None, exon=None,
        intron=None, cds_position=None, protein_position=None,
        max_af=None, popmax_af=None, nhomalt=None, bs1_threshold=None,
        revel=None, spliceai=None, clinvar_id=None, clinvar_sig=None,
        lof_mechanism=None, loeuf=None, pli=None,
        clinvar_missense_plp=None, clinvar_missense_blb=None, in_hotspot=None,
        human_pheno_score=None,
    )
    out = evaluate_all(ev)
    assert len(out) == 28
