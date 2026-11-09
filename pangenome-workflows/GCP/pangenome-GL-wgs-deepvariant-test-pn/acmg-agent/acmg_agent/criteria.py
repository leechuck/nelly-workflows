"""Deterministic evaluation of the 28 ACMG/AMP 2015 + ClinGen SVI criteria."""
from __future__ import annotations

from acmg_agent.inputs import LOF_CONSEQUENCES
from acmg_agent.model import Criterion, Evidence

REVEL_PP3 = 0.7
REVEL_BP4 = 0.2
SPLICEAI_PP3 = 0.5
SPLICEAI_BP4 = 0.1
BA1_AF = 0.05
BS1_DEFAULT = {"AD": 1e-4, "XD": 1e-4, "AR": 5e-3, "XR": 5e-3}
PM2_MAX_AF = {
    "AD": 1e-5,
    "XD": 1e-5,
    "AR": 1e-4,
    "XR": 1e-4,
    "XL": 1e-4,
    "MT": 1e-4,
    "OTHER": 1e-4,
}


def eval_pvs1(ev: Evidence) -> Criterion:
    csq = (ev.consequence or "").lower()
    if csq not in LOF_CONSEQUENCES and csq != "start_lost":
        return Criterion("PVS1", "not_met", None, f"consequence '{csq}' is not a canonical null variant")
    if ev.lof_mechanism is None:
        return Criterion("PVS1", "not_evaluable", None, "gene LoF mechanism unknown")
    if not ev.lof_mechanism:
        return Criterion("PVS1", "not_met", None, "LoF is not an established disease mechanism for this gene")
    if csq == "start_lost":
        return Criterion("PVS1", "met", "moderate", "start_lost in LoF-mechanism gene (PVS1_Moderate per Abou Tayoun 2018)")
    if ev.nmd_escape:
        return Criterion("PVS1", "met", "strong", f"{csq} in last exon / NMD-escape region; downgraded to PVS1_Strong")
    return Criterion("PVS1", "met", "very_strong", f"{csq} in LoF-mechanism gene (not in NMD-escape region)")


def eval_ps1(ev: Evidence) -> Criterion:
    if "missense" not in (ev.consequence or "").lower():
        return Criterion("PS1", "not_met", None, "not a missense variant")
    if ev.same_aa_ClinVar_PLP >= 2:
        return Criterion("PS1", "met", "strong", f"same amino-acid change in ClinVar P/LP with {ev.same_aa_ClinVar_PLP} stars")
    if ev.same_aa_ClinVar_PLP == 1:
        return Criterion("PS1", "met", "moderate", "same amino-acid change in ClinVar P/LP with 1 star")
    return Criterion("PS1", "not_met", None, "no qualifying same-AA ClinVar P/LP variant")


def eval_ps2(ev: Evidence) -> Criterion:
    return Criterion("PS2", "not_evaluable", None, "single-sample workflow: no parental genotypes for de novo verification")


def eval_ps3(ev: Evidence) -> Criterion:
    return Criterion("PS3", "not_evaluable", None, "no curated functional assay database wired")


def eval_ps4(ev: Evidence) -> Criterion:
    return Criterion("PS4", "not_evaluable", None, "no case-control odds-ratio table wired")


def eval_pm1(ev: Evidence) -> Criterion:
    if "missense" not in (ev.consequence or "").lower():
        return Criterion("PM1", "not_met", None, "not a missense variant")
    if ev.in_hotspot is None:
        return Criterion("PM1", "not_evaluable", None, "hotspot regions unavailable for gene")
    if ev.in_hotspot:
        return Criterion("PM1", "met", "moderate", "missense variant in curated mutational hotspot / critical domain")
    return Criterion("PM1", "not_met", None, "variant outside curated hotspot regions")


def eval_pm2(ev: Evidence) -> Criterion:
    thr = PM2_MAX_AF.get(ev.moi, 1e-4)
    if ev.max_af is None or ev.max_af == 0.0:
        return Criterion("PM2", "met", "supporting", f"absent from gnomAD (< {thr:.0e} for {ev.moi}; PM2_Supporting per ClinGen SVI)")
    if ev.max_af < thr:
        return Criterion("PM2", "met", "supporting", f"max AF {ev.max_af:.2e} < {thr:.0e} for {ev.moi} (PM2_Supporting per ClinGen SVI)")
    return Criterion("PM2", "not_met", None, f"max AF {ev.max_af:.2e} >= {thr:.0e} for {ev.moi}")


def eval_pm3(ev: Evidence) -> Criterion:
    if ev.moi not in ("AR", "XR"):
        return Criterion("PM3", "not_met", None, f"MOI is {ev.moi}, not recessive")
    if ev.zygosity == "hom":
        if ev.max_af is None or ev.max_af < 1e-3:
            return Criterion("PM3", "met", "moderate", "homozygous rare variant in recessive gene (1.0 pt -> PM3_Moderate)")
        return Criterion("PM3", "not_met", None, f"homozygous variant AF {ev.max_af:.2e} not rare enough for PM3")
    if ev.zygosity == "het":
        for v in ev.other_var_in_gene:
            if v.get("clinvar_sig") in ("Pathogenic", "Likely_pathogenic"):
                return Criterion("PM3", "met", "supporting", "unphased second P/LP variant in same recessive gene (0.5 pt -> PM3_Supporting)")
        return Criterion("PM3", "not_met", None, "no qualifying second P/LP variant in same recessive gene")
    return Criterion("PM3", "not_met", None, f"zygosity '{ev.zygosity}' does not satisfy PM3")


def eval_pm4(ev: Evidence) -> Criterion:
    csq = (ev.consequence or "").lower()
    if csq in ("inframe_insertion", "inframe_deletion"):
        if ev.in_repeat_region:
            return Criterion("PM4", "not_met", None, "in-frame indel in repeat region (covered by BP3)")
        return Criterion("PM4", "met", "moderate", f"{csq} outside repeat region")
    if csq == "stop_lost":
        return Criterion("PM4", "met", "moderate", "stop_lost altering protein length")
    return Criterion("PM4", "not_met", None, f"consequence '{csq}' is not an in-frame indel or stop_lost")


def eval_pm5(ev: Evidence) -> Criterion:
    if "missense" not in (ev.consequence or "").lower():
        return Criterion("PM5", "not_met", None, "not a missense variant")
    if ev.same_aa_ClinVar_PLP >= 1:
        return Criterion("PM5", "not_met", None, "same amino-acid change already covered by PS1")
    if ev.diff_aa_ClinVar_PLP >= 2:
        return Criterion("PM5", "met", "moderate", f"different pathogenic missense at same codon in ClinVar ({ev.diff_aa_ClinVar_PLP} stars)")
    if ev.diff_aa_ClinVar_PLP == 1:
        return Criterion("PM5", "met", "supporting", "different pathogenic missense at same codon in ClinVar (1 star)")
    return Criterion("PM5", "not_met", None, "no qualifying different-AA ClinVar P/LP at same codon")


def eval_pm6(ev: Evidence) -> Criterion:
    return Criterion("PM6", "not_evaluable", None, "single-sample workflow: no parental genotypes for assumed de novo check")


def eval_pp1(ev: Evidence) -> Criterion:
    return Criterion("PP1", "not_evaluable", None, "single-sample workflow: no family segregation data")


def eval_pp2(ev: Evidence) -> Criterion:
    if "missense" not in (ev.consequence or "").lower():
        return Criterion("PP2", "not_met", None, "not a missense variant")
    if ev.clinvar_missense_plp is None or ev.clinvar_missense_blb is None:
        return Criterion("PP2", "not_evaluable", None, "gene missense ClinVar counts unavailable")
    tot = ev.clinvar_missense_plp + ev.clinvar_missense_blb
    if ev.clinvar_missense_plp >= 20 and tot > 0 and (ev.clinvar_missense_plp / tot) >= 0.8:
        return Criterion("PP2", "met", "supporting", f"missense-constrained gene ({ev.clinvar_missense_plp}/{tot} ClinVar missense are P/LP)")
    return Criterion("PP2", "not_met", None, f"gene missense P/LP ratio ({ev.clinvar_missense_plp}/{tot}) below PP2 threshold")


def eval_pp3(ev: Evidence) -> Criterion:
    if eval_pvs1(ev).status == "met":
        return Criterion("PP3", "not_met", None, "not applied when PVS1 is met (avoids double-counting splice/LoF evidence)")
    if ev.revel is None and ev.spliceai is None:
        return Criterion("PP3", "not_evaluable", None, "no REVEL or SpliceAI score available")
    if (ev.revel is not None and ev.revel >= REVEL_PP3) or (ev.spliceai is not None and ev.spliceai >= SPLICEAI_PP3):
        return Criterion("PP3", "met", "supporting", f"in-silico damaging prediction (REVEL={ev.revel}, SpliceAI={ev.spliceai})")
    return Criterion("PP3", "not_met", None, f"in-silico below PP3 threshold (REVEL={ev.revel}, SpliceAI={ev.spliceai})")


def eval_pp4(ev: Evidence) -> Criterion:
    if ev.human_pheno_score is None:
        return Criterion("PP4", "not_evaluable", None, "no Exomiser human phenotype score available")
    if bool(ev.disease) and (
        ev.human_pheno_score >= 0.6
        or (ev.human_pheno_score >= 0.45 and bool(ev.pheno_evidence))
    ):
        return Criterion("PP4", "met", "supporting", f"specific phenotype match for {ev.disease} (score={ev.human_pheno_score:.2f})")
    return Criterion("PP4", "not_met", None, f"phenotype score ({ev.human_pheno_score:.2f}) or disease specificity below PP4 threshold")


def eval_pp5(ev: Evidence) -> Criterion:
    return Criterion("PP5", "not_applicable", None, "retired per ClinGen SVI (Biesecker et al. 2018)")


def eval_ba1(ev: Evidence) -> Criterion:
    if ev.max_af is None:
        return Criterion("BA1", "not_met", None, "variant absent from population databases")
    if ev.max_af >= BA1_AF:
        return Criterion("BA1", "met", "stand_alone", f"max AF {ev.max_af:.4f} >= {BA1_AF}")
    return Criterion("BA1", "not_met", None, f"max AF {ev.max_af:.2e} < {BA1_AF}")


def eval_bs1(ev: Evidence) -> Criterion:
    if ev.max_af is None:
        return Criterion("BS1", "not_met", None, "variant absent from population databases")
    if ev.max_af >= BA1_AF:
        return Criterion("BS1", "not_met", None, "superseded by BA1")
    thr = ev.bs1_threshold if ev.bs1_threshold is not None else BS1_DEFAULT.get(ev.moi, 5e-3)
    if ev.max_af >= thr:
        return Criterion("BS1", "met", "strong", f"max AF {ev.max_af:.4f} >= BS1 threshold {thr:.4f} for {ev.moi}")
    return Criterion("BS1", "not_met", None, f"max AF {ev.max_af:.2e} < BS1 threshold {thr:.4f}")


def eval_bs2(ev: Evidence) -> Criterion:
    if ev.nhomalt is not None:
        req = 1 if ev.moi in ("AD", "XD") else 2
        if ev.nhomalt >= req:
            return Criterion("BS2", "met", "strong", f"gnomAD nhomalt={ev.nhomalt} >= {req} for {ev.moi}")
        return Criterion("BS2", "not_met", None, f"gnomAD nhomalt={ev.nhomalt} < {req} for {ev.moi}")
    if ev.max_af is None or ev.max_af < 1e-4:
        return Criterion("BS2", "not_met", None, "variant absent or ultrarare in population databases")
    return Criterion("BS2", "not_evaluable", None, "gnomAD nhomalt unavailable (no slivar annotation)")


def eval_bs3(ev: Evidence) -> Criterion:
    return Criterion("BS3", "not_evaluable", None, "no curated benign functional assay database wired")


def eval_bs4(ev: Evidence) -> Criterion:
    return Criterion("BS4", "not_evaluable", None, "single-sample workflow: no family non-segregation data")


def eval_bp1(ev: Evidence) -> Criterion:
    if "missense" not in (ev.consequence or "").lower():
        return Criterion("BP1", "not_met", None, "not a missense variant")
    if ev.lof_mechanism is None or ev.clinvar_missense_plp is None:
        return Criterion("BP1", "not_evaluable", None, "gene mechanism or ClinVar counts unavailable")
    lof_plp = ev.clinvar_lof_plp if ev.clinvar_lof_plp is not None else 20
    tot = lof_plp + ev.clinvar_missense_plp
    if ev.lof_mechanism and lof_plp >= 10 and tot > 0 and (ev.clinvar_missense_plp / tot) <= 0.1:
        return Criterion("BP1", "met", "supporting", f"missense variant in truncating-only LoF gene (missense P/LP={ev.clinvar_missense_plp}, LoF P/LP={lof_plp})")
    return Criterion("BP1", "not_met", None, "gene is not exclusively truncating-variant pathogenic")


def eval_bp2(ev: Evidence) -> Criterion:
    if ev.moi in ("AD", "XD"):
        for v in ev.other_var_in_gene:
            if v.get("clinvar_sig") in ("Pathogenic", "Likely_pathogenic"):
                return Criterion("BP2", "met", "supporting", "co-occurrence with a known P/LP variant in the same dominant gene")
    return Criterion("BP2", "not_met", None, "no qualifying co-occurring P/LP variant in same dominant gene")


def eval_bp3(ev: Evidence) -> Criterion:
    csq = (ev.consequence or "").lower()
    if csq in ("inframe_insertion", "inframe_deletion") and ev.in_repeat_region:
        return Criterion("BP3", "met", "supporting", f"{csq} within a repetitive region without known function")
    return Criterion("BP3", "not_met", None, "not an in-frame indel in a repeat region")


def eval_bp4(ev: Evidence) -> Criterion:
    if ev.revel is None and ev.spliceai is None:
        return Criterion("BP4", "not_evaluable", None, "no REVEL or SpliceAI score available")
    csq = (ev.consequence or "").lower()
    if csq in LOF_CONSEQUENCES:
        return Criterion("BP4", "not_met", None, "canonical null variant is not eligible for BP4")
    if ev.revel is not None:
        if ev.revel <= REVEL_BP4 and (ev.spliceai is None or ev.spliceai < SPLICEAI_BP4):
            return Criterion("BP4", "met", "supporting", f"in-silico benign prediction (REVEL={ev.revel}, SpliceAI={ev.spliceai})")
        return Criterion("BP4", "not_met", None, f"in-silico above BP4 threshold (REVEL={ev.revel}, SpliceAI={ev.spliceai})")
    if ev.spliceai is not None and ev.spliceai < SPLICEAI_BP4 and csq in ("synonymous_variant", "intron_variant", "5_prime_utr_variant", "3_prime_utr_variant"):
        return Criterion("BP4", "met", "supporting", f"non-coding/synonymous variant with SpliceAI={ev.spliceai} < {SPLICEAI_BP4}")
    return Criterion("BP4", "not_met", None, f"in-silico scores do not establish BP4 (REVEL={ev.revel}, SpliceAI={ev.spliceai})")


def eval_bp5(ev: Evidence) -> Criterion:
    return Criterion("BP5", "not_met", None, "no alternate locus molecular cause established")


def eval_bp6(ev: Evidence) -> Criterion:
    return Criterion("BP6", "not_applicable", None, "retired per ClinGen SVI (Biesecker et al. 2018)")


def eval_bp7(ev: Evidence) -> Criterion:
    if (ev.consequence or "").lower() != "synonymous_variant":
        return Criterion("BP7", "not_met", None, "not a synonymous variant")
    if ev.spliceai is None:
        return Criterion("BP7", "not_evaluable", None, "synonymous variant without SpliceAI score")
    if ev.spliceai < SPLICEAI_BP4:
        return Criterion("BP7", "met", "supporting", f"synonymous variant with SpliceAI={ev.spliceai:.2f} < {SPLICEAI_BP4}")
    return Criterion("BP7", "not_met", None, f"synonymous variant with SpliceAI={ev.spliceai:.2f} >= {SPLICEAI_BP4}")


_EVALUATORS = (
    eval_pvs1, eval_ps1, eval_ps2, eval_ps3, eval_ps4,
    eval_pm1, eval_pm2, eval_pm3, eval_pm4, eval_pm5, eval_pm6,
    eval_pp1, eval_pp2, eval_pp3, eval_pp4, eval_pp5,
    eval_ba1, eval_bs1, eval_bs2, eval_bs3, eval_bs4,
    eval_bp1, eval_bp2, eval_bp3, eval_bp4, eval_bp5, eval_bp6, eval_bp7,
)


def evaluate_all(ev: Evidence) -> list[Criterion]:
    """Evaluate all 28 ACMG/AMP criteria in canonical order."""
    return [fn(ev) for fn in _EVALUATORS]
