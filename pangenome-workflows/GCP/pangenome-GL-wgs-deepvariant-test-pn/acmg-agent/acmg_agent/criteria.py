"""Deterministic evaluation of the 28 ACMG/AMP 2015 + ClinGen SVI criteria."""
from __future__ import annotations

from acmg_agent.inputs import LOF_CONSEQUENCES
from acmg_agent.model import Criterion, Evidence

# Pejaver et al. 2022 ClinGen SVI calibration for in-silico predictors:
# REVEL_PP3 uses >= thresholds; REVEL_BP4 uses <= thresholds.
REVEL_PP3: list[tuple[float, str]] = [
    (0.932, "strong"),
    (0.773, "moderate"),
    (0.644, "supporting"),
]
REVEL_BP4: list[tuple[float, str]] = [
    (0.003, "strong"),
    (0.016, "moderate"),
    (0.183, "supporting"),
]
SPLICEAI_PP3 = 0.2
SPLICEAI_BP4 = 0.1
ALPHAMISSENSE_PP3 = 0.564
ALPHAMISSENSE_BP4 = 0.34
LOEUF_PVS1_FALLBACK = 0.35

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

# ClinGen SVI PM3 point table (version 1.0) per proband allelic observation:
# - in trans with Pathogenic confirmed: 1.0 pt
# - phase unknown with Pathogenic: 0.5 pt
# - in trans with Likely pathogenic: 0.5 pt
# - phase unknown with Likely pathogenic: 0.25 pt
# - homozygous occurrence of rare variant: 0.5 pt
# Total points -> strength: 0.5..0.99 supporting, 1.0..1.99 moderate, 2.0..3.99 strong, >=4.0 very_strong
PM3_POINTS_TABLE = {
    "in_trans_pathogenic": 1.0,
    "unphased_pathogenic": 0.5,
    "in_trans_likely_pathogenic": 0.5,
    "unphased_likely_pathogenic": 0.25,
    "homozygous": 0.5,
}

_STRENGTH_ORDER = {
    "supporting": 1,
    "moderate": 2,
    "strong": 3,
    "very_strong": 4,
    "stand_alone": 5,
}


def _parse_exon_pair(exon: str | None) -> tuple[int, int] | None:
    if not exon or "/" not in exon:
        return None
    left, right = exon.split("/", 1)
    try:
        return (int(left), int(right))
    except ValueError:
        return None


def _removes_lt_10pct_in_last_exon(ev: Evidence) -> bool:
    pair = _parse_exon_pair(ev.exon)
    if pair is None:
        return False
    ex_idx, ex_cnt = pair
    if ex_cnt <= 1 or ex_idx != ex_cnt or ev.protein_position is None or ev.protein_position <= 0:
        return False
    if ev.protein_length is not None and ev.protein_length > 0:
        return ((ev.protein_length - ev.protein_position) / ev.protein_length) < 0.10
    # When total protein length is not separately annotated, a truncation in the terminal exon
    # of a >=10-exon transcript with a known protein position removes <10% of coding exons.
    return (1.0 / ex_cnt) <= 0.10


def eval_pvs1(ev: Evidence) -> Criterion:
    csq = (ev.consequence or "").lower()
    if csq not in LOF_CONSEQUENCES and csq != "start_lost":
        return Criterion("PVS1", "not_met", None, f"consequence '{csq}' is not a canonical null variant")

    lof_note = "LoF-mechanism gene"
    if ev.lof_mechanism is None:
        if ev.loeuf is not None and ev.loeuf < LOEUF_PVS1_FALLBACK:
            lof_note = f"lof_mechanism unknown; LOEUF={ev.loeuf:.2f} < {LOEUF_PVS1_FALLBACK} treated as LoF-mechanism gene"
        elif ev.loeuf is None:
            return Criterion("PVS1", "not_evaluable", None, "gene LoF mechanism and LOEUF unknown")
        else:
            return Criterion(
                "PVS1",
                "not_met",
                None,
                f"gene LoF mechanism unknown and LOEUF={ev.loeuf:.2f} >= {LOEUF_PVS1_FALLBACK}",
            )
    elif not ev.lof_mechanism:
        return Criterion("PVS1", "not_met", None, "LoF is not an established disease mechanism for this gene")

    if csq == "start_lost":
        return Criterion("PVS1", "met", "moderate", f"start_lost in {lof_note} (PVS1_Moderate per Abou Tayoun 2018)")
    if ev.nmd_escape or (_parse_exon_pair(ev.exon) is not None and _parse_exon_pair(ev.exon)[0] == _parse_exon_pair(ev.exon)[1] and _parse_exon_pair(ev.exon)[1] > 1):
        if _removes_lt_10pct_in_last_exon(ev):
            return Criterion(
                "PVS1",
                "met",
                "moderate",
                f"{csq} in last exon ({ev.exon}) at protein pos {ev.protein_position} removes <10% of protein in {lof_note}; PVS1_Moderate",
            )
        return Criterion("PVS1", "met", "strong", f"{csq} in last exon / NMD-escape region in {lof_note}; downgraded to PVS1_Strong")
    return Criterion("PVS1", "met", "very_strong", f"{csq} in {lof_note} (not in NMD-escape region)")


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


def _founder_af_context(ev: Evidence) -> str:
    if ev.max_af_any_pop is not None and (ev.max_af is None or ev.max_af_any_pop > ev.max_af):
        pop = ev.max_af_any_pop_pop or ev.max_af_any_pop_source or "founder"
        return f"; founder-population AF {ev.max_af_any_pop:.4f} ({pop})"
    return ""


def eval_pm2(ev: Evidence) -> Criterion:
    thr = PM2_MAX_AF.get(ev.moi, 1e-4)
    src = f" ({ev.max_af_source})" if ev.max_af_source else ""
    f_ctx = _founder_af_context(ev)
    if ev.max_af is None or ev.max_af == 0.0:
        return Criterion("PM2", "met", "supporting", f"absent from gnomAD (< {thr:.0e} for {ev.moi}; PM2_Supporting per ClinGen SVI{f_ctx})")
    if ev.max_af < thr:
        return Criterion("PM2", "met", "supporting", f"max AF {ev.max_af:.2e}{src} < {thr:.0e} for {ev.moi} (PM2_Supporting per ClinGen SVI{f_ctx})")
    return Criterion("PM2", "not_met", None, f"max AF {ev.max_af:.2e}{src} >= {thr:.0e} for {ev.moi}{f_ctx}")


def _pm3_strength_from_points(pts: float) -> str | None:
    if pts >= 4.0:
        return "very_strong"
    if pts >= 2.0:
        return "strong"
    if pts >= 1.0:
        return "moderate"
    if pts >= 0.5:
        return "supporting"
    return None


def eval_pm3(ev: Evidence) -> Criterion:
    if ev.moi not in ("AR", "XR"):
        return Criterion("PM3", "not_met", None, f"MOI is {ev.moi}, not recessive")
    if ev.zygosity == "hom":
        if ev.max_af is None or ev.max_af < 1e-3:
            pts = PM3_POINTS_TABLE["homozygous"]
            strength = _pm3_strength_from_points(pts)
            return Criterion(
                "PM3",
                "met",
                strength,
                f"homozygous rare variant in recessive gene ({pts:.1f} pts -> PM3_{strength.capitalize() if strength else 'None'} per ClinGen SVI)",
            )
        return Criterion("PM3", "not_met", None, f"homozygous variant AF {ev.max_af:.2e} not rare enough for PM3")
    if ev.zygosity == "het":
        best_pts = 0.0
        best_desc = ""
        for v in ev.other_var_in_gene:
            sig = v.get("clinvar_sig")
            in_trans = bool(v.get("phased_in_trans"))
            if sig == "Pathogenic":
                p = PM3_POINTS_TABLE["in_trans_pathogenic" if in_trans else "unphased_pathogenic"]
                desc = "in trans with Pathogenic" if in_trans else "phase unknown with Pathogenic"
            elif sig == "Likely_pathogenic":
                p = PM3_POINTS_TABLE["in_trans_likely_pathogenic" if in_trans else "unphased_likely_pathogenic"]
                desc = "in trans with Likely pathogenic" if in_trans else "phase unknown with Likely pathogenic"
            else:
                continue
            if p > best_pts:
                best_pts = p
                best_desc = desc
        strength = _pm3_strength_from_points(best_pts)
        if strength is not None:
            return Criterion(
                "PM3",
                "met",
                strength,
                f"second P/LP variant in same recessive gene: {best_desc} ({best_pts:.2g} pts -> PM3_{strength.capitalize()})",
            )
        if best_pts > 0:
            return Criterion("PM3", "not_met", None, f"second variant in recessive gene scored {best_pts:.2g} pts (< 0.5 pt threshold)")
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

    is_missense = "missense" in (ev.consequence or "").lower()
    am_val = ev.alphamissense if (ev.revel is None and is_missense) else None

    if ev.revel is None and ev.spliceai is None and am_val is None:
        return Criterion("PP3", "not_evaluable", None, "no REVEL, AlphaMissense, or SpliceAI score available")

    candidates: list[tuple[str, str]] = []
    if ev.revel is not None:
        for thr, strength in REVEL_PP3:
            if ev.revel >= thr:
                candidates.append((strength, f"REVEL={ev.revel:.3f} >= {thr} ({strength})"))
                break
    elif am_val is not None and am_val >= ALPHAMISSENSE_PP3:
        candidates.append(("supporting", f"AlphaMissense={am_val:.3f} >= {ALPHAMISSENSE_PP3} (supporting)"))

    if ev.spliceai is not None and ev.spliceai >= SPLICEAI_PP3:
        candidates.append(("supporting", f"SpliceAI={ev.spliceai:.2f} >= {SPLICEAI_PP3} (supporting)"))

    if candidates:
        best_strength, best_reason = max(candidates, key=lambda item: _STRENGTH_ORDER[item[0]])
        return Criterion("PP3", "met", best_strength, f"in-silico damaging prediction: {best_reason}")

    return Criterion(
        "PP3",
        "not_met",
        None,
        f"in-silico below PP3 thresholds (REVEL={ev.revel}, AlphaMissense={ev.alphamissense}, SpliceAI={ev.spliceai})",
    )


def eval_pp4(ev: Evidence) -> Criterion:
    if ev.human_pheno_score is None:
        return Criterion("PP4", "not_evaluable", None, "no Exomiser human phenotype score available")
    if bool(ev.disease) and ev.human_pheno_score >= 0.6:
        return Criterion("PP4", "met", "supporting", f"specific phenotype match for {ev.disease} (score={ev.human_pheno_score:.2f} >= 0.60)")
    return Criterion("PP4", "not_met", None, f"phenotype score ({ev.human_pheno_score:.2f}) or disease specificity below PP4 threshold (0.60)")


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
    csq = (ev.consequence or "").lower()
    is_missense = "missense" in csq
    am_val = ev.alphamissense if (ev.revel is None and is_missense) else None

    if ev.revel is None and ev.spliceai is None and am_val is None:
        return Criterion("BP4", "not_evaluable", None, "no REVEL, AlphaMissense, or SpliceAI score available")
    if csq in LOF_CONSEQUENCES:
        return Criterion("BP4", "not_met", None, "canonical null variant is not eligible for BP4")

    # Any present predictor exceeding the BP4 benign cutoff prevents BP4 from being met
    if ev.spliceai is not None and ev.spliceai > SPLICEAI_BP4:
        return Criterion("BP4", "not_met", None, f"SpliceAI={ev.spliceai:.2f} > {SPLICEAI_BP4}")

    candidates: list[tuple[str, str]] = []
    if ev.revel is not None:
        matched_revel = False
        for thr, strength in REVEL_BP4:
            if ev.revel <= thr:
                candidates.append((strength, f"REVEL={ev.revel:.3f} <= {thr} ({strength})"))
                matched_revel = True
                break
        if not matched_revel:
            return Criterion("BP4", "not_met", None, f"REVEL={ev.revel:.3f} > {REVEL_BP4[-1][0]}")
    elif am_val is not None:
        if am_val <= ALPHAMISSENSE_BP4:
            candidates.append(("supporting", f"AlphaMissense={am_val:.3f} <= {ALPHAMISSENSE_BP4} (supporting)"))
        else:
            return Criterion("BP4", "not_met", None, f"AlphaMissense={am_val:.3f} > {ALPHAMISSENSE_BP4}")

    if ev.spliceai is not None and ev.spliceai <= SPLICEAI_BP4:
        if ev.revel is not None or am_val is not None or csq in ("synonymous_variant", "intron_variant", "5_prime_utr_variant", "3_prime_utr_variant"):
            candidates.append(("supporting", f"SpliceAI={ev.spliceai:.2f} <= {SPLICEAI_BP4} (supporting)"))

    if candidates:
        best_strength, best_reason = max(candidates, key=lambda item: _STRENGTH_ORDER[item[0]])
        return Criterion("BP4", "met", best_strength, f"in-silico benign prediction: {best_reason}")

    return Criterion(
        "BP4",
        "not_met",
        None,
        f"in-silico scores do not establish BP4 (REVEL={ev.revel}, AlphaMissense={ev.alphamissense}, SpliceAI={ev.spliceai})",
    )


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
