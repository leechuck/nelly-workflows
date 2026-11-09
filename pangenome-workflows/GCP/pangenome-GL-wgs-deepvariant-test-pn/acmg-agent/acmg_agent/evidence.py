"""Evidence builder combining Exomiser candidates, VEP CSQ, ClinVar, gene table, and slivar."""
from __future__ import annotations

from acmg_agent.inputs import (
    ALL_POP_AF_FIELDS,
    ClinvarIndex,
    NON_BOTTLENECKED_AF_FIELDS,
    VepAnnotation,
    VepIndex,
    extract_pop_code,
    norm_chrom,
    parse_exomiser_freq,
    parse_float,
    parse_int,
    _normalize_clnsig,
)
from acmg_agent.model import Evidence


def _resolve_population_frequencies(
    vep_hit: VepAnnotation | None,
    slivar_popmax: float | None,
    exomiser_af: float | None,
) -> tuple[float | None, str | None, float | None, str | None, str | None]:
    """Return (max_af, max_af_source, max_af_any_pop, max_af_any_pop_source, max_af_any_pop_pop).

    ClinGen SVI PM2/BS1/BA1 calibration uses the maximum AF across non-bottlenecked reference
    populations (AFR, AMR, EAS, NFE, SAS across gnomADe_* and gnomADg_*), excluding bottlenecked
    or founder populations (ASJ, FIN, AMI, MID) and REMAINING/OTH buckets.
    """
    max_af: float | None = None
    max_af_source: str | None = None

    if vep_hit is not None:
        non_bottle = [
            (k, vep_hit.pop_afs[k])
            for k in NON_BOTTLENECKED_AF_FIELDS
            if k in vep_hit.pop_afs and vep_hit.pop_afs[k] is not None
        ]
        if non_bottle:
            max_af_source, max_af = max(non_bottle, key=lambda item: item[1])
        else:
            overall = [
                (k, v)
                for k, v in (("gnomADe_AF", vep_hit.gnomade_af), ("gnomADg_AF", vep_hit.gnomadg_af))
                if v is not None
            ]
            if overall:
                max_af_source, max_af = max(overall, key=lambda item: item[1])

    if max_af is None and slivar_popmax is not None:
        max_af = slivar_popmax
        max_af_source = "gnomad_popmax_af"

    if max_af is None and vep_hit is not None and vep_hit.max_af is not None:
        max_af = vep_hit.max_af
        max_af_source = "MAX_AF"

    if max_af is None and exomiser_af is not None:
        max_af = exomiser_af
        max_af_source = "MAX_FREQ"

    # Founder / all-population maximum for clinical context
    any_af: float | None = None
    any_src: str | None = None
    any_pop: str | None = None
    if vep_hit is not None:
        all_pops = [
            (k, vep_hit.pop_afs[k])
            for k in ALL_POP_AF_FIELDS
            if k in vep_hit.pop_afs and vep_hit.pop_afs[k] is not None
        ]
        if all_pops:
            any_src, any_af = max(all_pops, key=lambda item: item[1])
            any_pop = extract_pop_code(any_src)
        if vep_hit.max_af is not None and (any_af is None or vep_hit.max_af > any_af):
            any_af = vep_hit.max_af
            any_src = vep_hit.max_af_pops or "MAX_AF"
            if vep_hit.max_af_pops:
                any_pop = extract_pop_code(vep_hit.max_af_pops)

    if any_af is None:
        any_af = max_af
        any_src = max_af_source
        if max_af_source in ALL_POP_AF_FIELDS:
            any_pop = extract_pop_code(max_af_source)

    return (max_af, max_af_source, any_af, any_src, any_pop)


def _infer_zygosity(gt: str, chrom: str, moi: str) -> str:
    g = (gt or "").strip()
    if not g or g in (".", "./.", ".|."):
        return "unknown"
    c = norm_chrom(chrom).upper()
    if g in ("1", "0|1:hemi"):
        return "hemi"
    sep = "|" if "|" in g else "/"
    parts = [p for p in g.split(sep) if p != ""]
    if len(parts) == 1 and parts[0] not in ("0", "."):
        return "hemi"
    if len(parts) == 2:
        a, b = parts
        if a == b and a not in ("0", "."):
            if c in ("X", "Y") and moi in ("XL", "XR", "XD"):
                return "hemi"
            return "hom"
        if a != b and ("0" in (a, b) or (a not in ("0", ".") and b not in ("0", "."))):
            return "het"
    return "unknown"


def _parse_hgvs_parts(hgvs: str) -> tuple[str | None, str | None, str | None]:
    """Parse Exomiser HGVS field like 'HEXA:NM_000520.6:c.1421+1G>C:p.?' into (transcript, hgvsc, hgvsp)."""
    if not hgvs:
        return (None, None, None)
    parts = [p.strip() for p in hgvs.split(":")]
    tx = None
    hgvsc = None
    hgvsp = None
    for p in parts:
        if p.startswith(("NM_", "NR_", "ENST")):
            tx = p
        elif p.startswith(("c.", "n.", "m.")):
            hgvsc = p
        elif p.startswith("p."):
            hgvsp = p
    return (tx, hgvsc, hgvsp)


def _check_in_hotspot(chrom: str, pos: int, regions: list[tuple[str, int, int]]) -> bool:
    if not regions:
        return False
    c = norm_chrom(chrom)
    for r_chrom, start, end in regions:
        if r_chrom == c and start <= pos <= end:
            return True
    return False


def _collect_other_vars_in_gene(
    gene: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    variants_rows: list[dict],
    clinvar: ClinvarIndex | None,
) -> list[dict]:
    if not gene or not variants_rows:
        return []
    c_norm = norm_chrom(chrom)
    other: list[dict] = []
    seen: set[tuple[str, int, str, str]] = set()
    for r in variants_rows:
        r_gene = (r.get("GENE_SYMBOL") or r.get("gene") or "").strip()
        if r_gene != gene:
            continue
        r_chrom = norm_chrom(r.get("CONTIG") or r.get("chrom") or "")
        r_pos = parse_int(r.get("START") or r.get("pos")) or 0
        r_ref = (r.get("REF") or r.get("ref") or "").strip()
        r_alt = (r.get("ALT") or r.get("alt") or "").strip()
        if r_chrom == c_norm and r_pos == pos and r_ref == ref and r_alt == alt:
            continue
        key = (r_chrom, r_pos, r_ref, r_alt)
        if key in seen:
            continue
        seen.add(key)

        cv_sig = None
        cv_stars = 0
        if clinvar is not None and r_chrom and r_pos:
            hit = clinvar.same_variant(r_chrom, r_pos, r_ref, r_alt)
            if hit is not None:
                cv_sig = hit.clinvar_sig
                cv_stars = hit.clinvar_stars
        if cv_sig is None:
            raw_cv = r.get("CLINVAR_PRIMARY_INTERPRETATION") or r.get("clinvar_sig") or r.get("clinvar") or ""
            cv_sig = _normalize_clnsig(str(raw_cv))

        other.append({
            "chrom": r_chrom,
            "pos": r_pos,
            "ref": r_ref,
            "alt": r_alt,
            "gt": (r.get("GENOTYPE") or r.get("gt") or "").strip(),
            "consequence": (r.get("FUNCTIONAL_CLASS") or r.get("consequence") or "").strip().lower(),
            "clinvar_sig": cv_sig,
            "clinvar_stars": cv_stars,
            "max_af": parse_exomiser_freq(r.get("MAX_FREQ")) if "MAX_FREQ" in r else parse_float(r.get("max_af")),
        })
    return other


def build_evidence(
    candidate: dict,
    variants_rows: list[dict],
    vep: VepIndex | None,
    clinvar: ClinvarIndex | None,
    gene_table: dict[str, dict] | None,
    slivar_info: dict[tuple[str, int, str, str], dict] | None = None,
) -> Evidence:
    """Assemble an Evidence record for a candidate variant from all deterministic sources."""
    gene = (candidate.get("gene") or "").strip()
    chrom = str(candidate.get("chrom") or "").strip()
    pos = parse_int(candidate.get("pos")) or 0
    ref = str(candidate.get("ref") or "").strip()
    alt = str(candidate.get("alt") or "").strip()

    gt_row = gene_table.get(gene) if gene_table and gene else None
    moi = (candidate.get("moi") or (gt_row["moi"] if gt_row and gt_row.get("moi") else None) or "OTHER").strip()

    vep_hit = vep.lookup(chrom, pos, ref, alt, gene=gene or None) if vep is not None and chrom and pos else None
    raw_gt = (candidate.get("gt") or (vep_hit.gt if vep_hit else "") or "").strip()
    zygosity = _infer_zygosity(raw_gt, chrom, moi)

    tx_fallback, hgvsc_fallback, hgvsp_fallback = _parse_hgvs_parts(str(candidate.get("hgvs") or ""))
    transcript = (vep_hit.transcript if vep_hit and vep_hit.transcript else tx_fallback)
    mane_select = bool(vep_hit.mane_select) if vep_hit else False
    consequence = (
        (vep_hit.consequence if vep_hit and vep_hit.consequence else str(candidate.get("effect") or "").strip().lower())
        or "unknown"
    )
    hgvsc = (vep_hit.hgvsc if vep_hit and vep_hit.hgvsc else hgvsc_fallback)
    hgvsp = (vep_hit.hgvsp if vep_hit and vep_hit.hgvsp else hgvsp_fallback)
    exon = vep_hit.exon if vep_hit else None
    intron = vep_hit.intron if vep_hit else None
    cds_position = vep_hit.cds_position if vep_hit else None
    protein_position = vep_hit.protein_position if vep_hit else None
    protein_length = vep_hit.protein_length if vep_hit else None
    nmd_escape = bool(vep_hit.nmd_escape) if vep_hit else False
    in_repeat_region = bool(vep_hit.in_repeat_region) if vep_hit else False

    revel = vep_hit.revel if vep_hit else None
    spliceai = vep_hit.spliceai if vep_hit else None
    alphamissense = vep_hit.alphamissense if vep_hit else None

    # Population frequency: Exomiser MAX_FREQ is a percentage (0..100), converted once via parse_exomiser_freq
    exomiser_af = parse_exomiser_freq(candidate.get("max_freq"))
    slivar_key = (norm_chrom(chrom), pos, ref, alt)
    slivar_rec = slivar_info.get(slivar_key) if slivar_info else None
    slivar_popmax = slivar_rec.get("popmax_af") if slivar_rec else None
    nhomalt = slivar_rec.get("nhomalt") if slivar_rec else None

    (
        max_af,
        max_af_source,
        max_af_any_pop,
        max_af_any_pop_source,
        max_af_any_pop_pop,
    ) = _resolve_population_frequencies(vep_hit, slivar_popmax, exomiser_af)
    popmax_af = slivar_popmax if slivar_popmax is not None else (vep_hit.max_af if vep_hit else None)

    # ClinVar
    cv_hit = clinvar.same_variant(chrom, pos, ref, alt) if clinvar is not None and chrom and pos else None
    clinvar_id = cv_hit.clinvar_id if cv_hit else None
    clinvar_sig = cv_hit.clinvar_sig if cv_hit else None
    clinvar_stars = cv_hit.clinvar_stars if cv_hit else 0

    same_aa_plp = 0
    diff_aa_plp = 0
    if clinvar is not None and chrom and pos and hgvsp:
        same_aa_plp, diff_aa_plp = clinvar.codon_plp_stars(chrom, pos, ref, alt, hgvsp)

    # Gene-level
    if gt_row is not None:
        lof_mechanism = gt_row["lof_mechanism"]
        loeuf = gt_row["loeuf"]
        pli = gt_row["pli"]
        clinvar_missense_plp = gt_row["clinvar_missense_plp"]
        clinvar_missense_blb = gt_row["clinvar_missense_blb"]
        clinvar_lof_plp = gt_row["clinvar_lof_plp"]
        bs1_threshold = gt_row["bs1_threshold"]
        in_hotspot = _check_in_hotspot(chrom, pos, gt_row["hotspot_regions"])
    else:
        lof_mechanism = None
        loeuf = None
        pli = None
        clinvar_missense_plp = None
        clinvar_missense_blb = None
        clinvar_lof_plp = None
        bs1_threshold = None
        in_hotspot = None

    other_vars = _collect_other_vars_in_gene(gene, chrom, pos, ref, alt, variants_rows, clinvar)
    human_pheno_score = parse_float(candidate.get("pheno_score"))
    if human_pheno_score is None:
        human_pheno_score = parse_float(candidate.get("human_pheno_score"))

    return Evidence(
        gene=gene,
        moi=moi,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        zygosity=zygosity,
        transcript=transcript,
        mane_select=mane_select,
        consequence=consequence,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        exon=exon,
        intron=intron,
        cds_position=cds_position,
        protein_position=protein_position,
        nmd_escape=nmd_escape,
        in_repeat_region=in_repeat_region,
        max_af=max_af,
        popmax_af=popmax_af,
        nhomalt=nhomalt,
        bs1_threshold=bs1_threshold,
        revel=revel,
        spliceai=spliceai,
        clinvar_id=clinvar_id,
        clinvar_sig=clinvar_sig,
        clinvar_stars=clinvar_stars,
        same_aa_ClinVar_PLP=same_aa_plp,
        diff_aa_ClinVar_PLP=diff_aa_plp,
        lof_mechanism=lof_mechanism,
        loeuf=loeuf,
        pli=pli,
        clinvar_missense_plp=clinvar_missense_plp,
        clinvar_missense_blb=clinvar_missense_blb,
        in_hotspot=in_hotspot,
        other_var_in_gene=other_vars,
        human_pheno_score=human_pheno_score,
        pheno_evidence=str(candidate.get("pheno_evidence") or ""),
        disease=str(candidate.get("disease") or ""),
        clinvar_lof_plp=clinvar_lof_plp,
        alphamissense=alphamissense,
        protein_length=protein_length,
        max_af_source=max_af_source,
        max_af_any_pop=max_af_any_pop,
        max_af_any_pop_source=max_af_any_pop_source,
        max_af_any_pop_pop=max_af_any_pop_pop,
    )
