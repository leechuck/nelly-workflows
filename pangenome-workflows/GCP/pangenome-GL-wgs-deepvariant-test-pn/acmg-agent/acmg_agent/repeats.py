"""STRchive-backed short tandem repeat classification for ExpansionHunter calls."""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from acmg_agent.inputs import norm_chrom, parse_int


@dataclass
class Locus:
    gene: str
    chrom: str
    start_hg38: int
    stop_hg38: int
    motif: str
    inheritance: str
    disease: str
    benign_min: int | None
    benign_max: int | None
    intermediate_min: int | None
    intermediate_max: int | None
    pathogenic_min: int | None
    pathogenic_max: int | None
    all_regions: list[tuple[str, int, int]] = field(default_factory=list)


@dataclass
class RepeatCall:
    gene: str
    genotype: str
    max_allele: int
    tier: str                     # "normal" | "intermediate" | "pathogenic" | "unknown_locus"
    ranges: dict
    disease: str
    moi: str
    ci: str


def load_strchive(path: Path | str) -> dict[str, Locus]:
    """Load STRchive locus definitions keyed by gene symbol."""
    p = Path(path)
    raw_list = json.loads(p.read_text())
    loci: dict[str, Locus] = {}
    for entry in raw_list:
        gene = str(entry.get("gene") or "").strip()
        if not gene:
            continue
        chrom = str(entry.get("chrom") or "").strip()
        start_hg38 = parse_int(entry.get("start_hg38")) or 0
        stop_hg38 = parse_int(entry.get("stop_hg38")) or 0
        motifs = (
            entry.get("pathogenic_motif_reference_orientation")
            or entry.get("reference_motif_reference_orientation")
            or []
        )
        motif = str(motifs[0]) if motifs else ""
        inh_raw = entry.get("inheritance") or []
        inheritance = "/".join(str(x) for x in inh_raw) if isinstance(inh_raw, list) else str(inh_raw)
        disease = str(entry.get("disease") or "").strip()

        loc = Locus(
            gene=gene,
            chrom=chrom,
            start_hg38=start_hg38,
            stop_hg38=stop_hg38,
            motif=motif,
            inheritance=inheritance,
            disease=disease,
            benign_min=parse_int(entry.get("benign_min")),
            benign_max=parse_int(entry.get("benign_max")),
            intermediate_min=parse_int(entry.get("intermediate_min")),
            intermediate_max=parse_int(entry.get("intermediate_max")),
            pathogenic_min=parse_int(entry.get("pathogenic_min")),
            pathogenic_max=parse_int(entry.get("pathogenic_max")),
            all_regions=[(chrom, start_hg38, stop_hg38)],
        )
        if gene in loci:
            loci[gene].all_regions.append((chrom, start_hg38, stop_hg38))
        else:
            loci[gene] = loc
    return loci


def _parse_region(region: str) -> tuple[str, int, int] | None:
    if not region or ":" not in region or "-" not in region:
        return None
    chrom_part, span = region.strip().split(":", 1)
    start_s, stop_s = span.split("-", 1)
    try:
        return (norm_chrom(chrom_part), int(start_s), int(stop_s))
    except ValueError:
        return None


def _overlaps_locus(region_str: str, locus: Locus) -> bool:
    parsed = _parse_region(region_str)
    if parsed is None:
        return True
    r_chrom, r_start, r_stop = parsed
    regions = locus.all_regions or [(locus.chrom, locus.start_hg38, locus.stop_hg38)]
    for l_chrom, l_start, l_stop in regions:
        if norm_chrom(l_chrom) == r_chrom and max(r_start, l_start) <= min(r_stop, l_stop):
            return True
    return False


def classify_repeat(eh_variant: dict, loci: dict[str, Locus]) -> RepeatCall:
    """Classify an ExpansionHunter repeat variant against STRchive locus definitions."""
    raw_gene = str(
        eh_variant.get("VariantId")
        or eh_variant.get("LocusId")
        or eh_variant.get("gene")
        or ""
    ).strip()
    gt = str(
        eh_variant.get("Genotype")
        or eh_variant.get("genotype")
        or eh_variant.get("genotypes")
        or ""
    ).strip()
    ci = str(
        eh_variant.get("GenotypeConfidenceInterval")
        or eh_variant.get("ci")
        or ""
    ).strip()
    region = str(
        eh_variant.get("ReferenceRegion")
        or eh_variant.get("region")
        or ""
    ).strip()

    alleles = [int(a) for a in gt.split("/") if a.strip().lstrip("-").isdigit()]
    max_allele = max(alleles) if alleles else 0
    min_allele = min(alleles) if alleles else 0

    locus = loci.get(raw_gene)
    if locus is None and raw_gene:
        upper_map = {k.upper(): v for k, v in loci.items()}
        locus = upper_map.get(raw_gene.upper())

    if locus is None or (region and not _overlaps_locus(region, locus)):
        return RepeatCall(
            gene=raw_gene,
            genotype=gt,
            max_allele=max_allele,
            tier="unknown_locus",
            ranges={},
            disease="",
            moi="",
            ci=ci,
        )

    ranges = {
        "benign": (locus.benign_min, locus.benign_max),
        "intermediate": (locus.intermediate_min, locus.intermediate_max),
        "pathogenic": (locus.pathogenic_min, locus.pathogenic_max),
    }

    # For autosomal recessive loci (e.g. FXN), full biallelic pathogenicity uses min_allele
    # when two alleles are called, whereas intermediate/carrier or dominant loci use max_allele.
    eff_path_allele = min_allele if (locus.inheritance == "AR" and len(alleles) == 2) else max_allele

    if locus.pathogenic_min is not None and eff_path_allele >= locus.pathogenic_min:
        tier = "pathogenic"
    elif locus.intermediate_min is not None and max_allele >= locus.intermediate_min:
        tier = "intermediate"
    elif locus.benign_max is not None and max_allele > locus.benign_max:
        tier = "intermediate"
    else:
        tier = "normal"

    return RepeatCall(
        gene=locus.gene.upper() if raw_gene.isupper() else locus.gene,
        genotype=gt,
        max_allele=max_allele,
        tier=tier,
        ranges=ranges,
        disease=locus.disease,
        moi=locus.inheritance,
        ci=ci,
    )
