"""Input readers and indexes for VEP CSQ VCFs, ClinVar tabix VCFs, gene tables, and slivar VCFs."""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
import pysam


LOF_CONSEQUENCES = frozenset({
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
})

_CODON_RE = re.compile(r"p\.([A-Za-z]{1,3})(\d+)([A-Za-z=*?]{1,3})")


def norm_chrom(chrom: str) -> str:
    """Return bare chromosome name without leading 'chr' (except chrM -> MT/M preserved via helpers)."""
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        return c[3:]
    return c


def parse_float(val: object) -> float | None:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in (".", "NA", "nan", "None"):
        return None
    # Handle multi-valued fields separated by '&', ',', or '|'
    for sep in ("&", ",", "|"):
        if sep in s:
            parts = [parse_float(p) for p in s.split(sep)]
            valid = [p for p in parts if p is not None]
            return max(valid) if valid else None
    try:
        return float(s)
    except ValueError:
        return None


def parse_int(val: object) -> int | None:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in (".", "NA", "nan", "None"):
        return None
    # Handle ranges such as "1497-1514" or "500-505"
    if "-" in s:
        s = s.split("-", 1)[0].strip()
    if "/" in s:
        s = s.split("/", 1)[0].strip()
    try:
        return int(s)
    except ValueError:
        return None


def parse_exomiser_freq(val: object) -> float | None:
    """Convert Exomiser MAX_FREQ (reported as percentage 0..100) into a fraction (0..1)."""
    f = parse_float(val)
    if f is None:
        return None
    return f / 100.0


@dataclass
class VepAnnotation:
    gene: str
    transcript: str | None
    mane_select: bool
    canonical: bool
    consequence: str
    hgvsc: str | None
    hgvsp: str | None
    exon: str | None
    intron: str | None
    cds_position: int | None
    protein_position: int | None
    nmd_escape: bool
    in_repeat_region: bool
    max_af: float | None
    gnomade_af: float | None
    gnomadg_af: float | None
    revel: float | None
    spliceai: float | None
    gt: str | None = None


def _select_consequence(raw_csq: str) -> str:
    if not raw_csq:
        return "unknown"
    parts = [p.strip().lower() for p in raw_csq.split("&") if p.strip()]
    for p in parts:
        if p in LOF_CONSEQUENCES:
            return p
    return parts[0] if parts else "unknown"


def _compute_nmd_escape(consequence: str, exon: str | None) -> bool:
    if consequence not in ("stop_gained", "frameshift_variant"):
        return False
    if not exon or "/" not in exon:
        return False
    left, right = exon.split("/", 1)
    try:
        ex_num = int(left)
        ex_tot = int(right)
        return ex_tot > 1 and ex_num == ex_tot
    except ValueError:
        return False


def _extract_spliceai(entry: dict[str, str]) -> float | None:
    vals: list[float] = []
    for k in ("SpliceAI", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"):
        if k in entry and entry[k]:
            f = parse_float(entry[k])
            if f is not None:
                vals.append(f)
    return max(vals) if vals else None


class VepIndex:
    """Index of VEP CSQ annotations keyed by normalized (chrom, pos, ref, alt)."""

    def __init__(self, records: dict[tuple[str, int, str, str], list[VepAnnotation]]):
        self._records = records

    @classmethod
    def from_vcf(cls, path: Path | str) -> "VepIndex":
        vcf_path = Path(path)
        records: dict[tuple[str, int, str, str], list[VepAnnotation]] = {}
        with pysam.VariantFile(str(vcf_path)) as vf:
            csq_fields: list[str] = []
            if "CSQ" in vf.header.info:
                desc = vf.header.info["CSQ"].description or ""
                if "Format:" in desc:
                    fmt = desc.split("Format:", 1)[1].strip().rstrip('"')
                    csq_fields = [f.strip() for f in fmt.split("|")]

            for rec in vf:
                chrom = norm_chrom(rec.chrom)
                pos = int(rec.pos)
                ref = str(rec.ref)
                alts = rec.alts or ()
                gt_str: str | None = None
                if len(rec.samples) > 0:
                    s0 = rec.samples[0]
                    gt_tuple = s0.get("GT")
                    if gt_tuple is not None:
                        sep = "|" if s0.phased else "/"
                        gt_str = sep.join("." if a is None else str(a) for a in gt_tuple)

                csq_raw = rec.info.get("CSQ", ())
                if isinstance(csq_raw, str):
                    csq_raw = (csq_raw,)

                parsed_annots: list[VepAnnotation] = []
                for item in csq_raw:
                    vals = str(item).split("|")
                    entry = dict(zip(csq_fields, vals)) if csq_fields else {}
                    consequence = _select_consequence(entry.get("Consequence", ""))
                    exon = entry.get("EXON") or None
                    intron = entry.get("INTRON") or None
                    mane = bool(entry.get("MANE_SELECT", "").strip())
                    canon = entry.get("CANONICAL", "").strip().upper() == "YES"
                    repeat_flag = bool(entry.get("Repeat", "").strip()) if "Repeat" in entry else False
                    annot = VepAnnotation(
                        gene=entry.get("SYMBOL", "").strip(),
                        transcript=entry.get("Feature", "").strip() or None,
                        mane_select=mane,
                        canonical=canon,
                        consequence=consequence,
                        hgvsc=entry.get("HGVSc", "").strip() or None,
                        hgvsp=entry.get("HGVSp", "").strip() or None,
                        exon=exon,
                        intron=intron,
                        cds_position=parse_int(entry.get("CDS_position")),
                        protein_position=parse_int(entry.get("Protein_position")),
                        nmd_escape=_compute_nmd_escape(consequence, exon),
                        in_repeat_region=repeat_flag,
                        max_af=parse_float(entry.get("MAX_AF")),
                        gnomade_af=parse_float(entry.get("gnomADe_AF")),
                        gnomadg_af=parse_float(entry.get("gnomADg_AF")),
                        revel=parse_float(entry.get("REVEL")),
                        spliceai=_extract_spliceai(entry),
                        gt=gt_str,
                    )
                    parsed_annots.append(annot)

                for alt in alts:
                    key = (chrom, pos, ref, str(alt))
                    records.setdefault(key, []).extend(parsed_annots)

        return cls(records)

    def lookup(self, chrom: str, pos: int, ref: str, alt: str, gene: str | None = None) -> VepAnnotation | None:
        key = (norm_chrom(chrom), int(pos), str(ref), str(alt))
        annots = self._records.get(key)
        if not annots:
            return None
        pool = [a for a in annots if gene and a.gene == gene] or annots
        for a in pool:
            if a.mane_select:
                return a
        for a in pool:
            if a.canonical:
                return a
        return pool[0]


@dataclass
class ClinvarHit:
    clinvar_id: str | None
    clinvar_sig: str | None
    clinvar_stars: int


def _normalize_clnsig(raw: str) -> str | None:
    if not raw or raw == ".":
        return None
    s = raw.strip()
    low = s.lower()
    if "conflicting" in low:
        return "Conflicting"
    # Handle slash-separated combined significance like Pathogenic/Likely_pathogenic
    if "pathogenic" in low and "likely_pathogenic" in low:
        return "Pathogenic"
    if low.startswith("pathogenic"):
        return "Pathogenic"
    if low.startswith("likely_pathogenic"):
        return "Likely_pathogenic"
    if "benign" in low and "likely_benign" in low:
        return "Benign"
    if low.startswith("benign"):
        return "Benign"
    if low.startswith("likely_benign"):
        return "Likely_benign"
    if "uncertain_significance" in low:
        return "Uncertain_significance"
    return s


def _clnrevstat_to_stars(raw: str) -> int:
    if not raw:
        return 0
    low = raw.strip().lower()
    if "practice_guideline" in low:
        return 4
    if "reviewed_by_expert_panel" in low:
        return 3
    if "criteria_provided" in low and "_multiple_submitters" in low and "_no_conflicts" in low:
        return 2
    if "criteria_provided" in low and "_single_submitter" in low:
        return 1
    return 0


def _parse_vcf_info(info_str: str) -> dict[str, str]:
    res: dict[str, str] = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            res[k] = v
        elif part:
            res[part] = "True"
    return res


def _parse_codon(hgvsp: str | None) -> tuple[str, int, str] | None:
    if not hgvsp:
        return None
    m = _CODON_RE.search(hgvsp)
    if not m:
        return None
    ref_aa, pos_s, alt_aa = m.group(1), m.group(2), m.group(3)
    if alt_aa in ("?", "=", "*", "Ter", "fs"):
        return None
    return (ref_aa, int(pos_s), alt_aa)


class ClinvarIndex:
    """Tabix-backed ClinVar lookup for exact variant and same-codon P/LP alleles."""

    def __init__(self, path: Path | str):
        self.path = Path(path)
        self.tbx = pysam.TabixFile(str(self.path))
        self._contigs = set(self.tbx.contigs)

    def _resolve_contig(self, chrom: str) -> str | None:
        c = str(chrom).strip()
        if c in self._contigs:
            return c
        bare = norm_chrom(c)
        if bare in self._contigs:
            return bare
        with_chr = f"chr{bare}"
        if with_chr in self._contigs:
            return with_chr
        return None

    def same_variant(self, chrom: str, pos: int, ref: str, alt: str) -> ClinvarHit | None:
        contig = self._resolve_contig(chrom)
        if contig is None:
            return None
        p = int(pos)
        try:
            rows = list(self.tbx.fetch(contig, max(0, p - 1), p))
        except ValueError:
            return None
        for line in rows:
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            r_pos = int(cols[1])
            r_id = cols[2] if cols[2] != "." else None
            r_ref = cols[3]
            r_alts = cols[4].split(",")
            if r_pos == p and r_ref == ref and alt in r_alts:
                info = _parse_vcf_info(cols[7])
                sig = _normalize_clnsig(info.get("CLNSIG", ""))
                stars = _clnrevstat_to_stars(info.get("CLNREVSTAT", ""))
                return ClinvarHit(clinvar_id=r_id, clinvar_sig=sig, clinvar_stars=stars)
        return None

    def codon_plp_stars(self, chrom: str, pos: int, ref: str, alt: str, hgvsp: str | None) -> tuple[int, int]:
        """Return (same_aa_max_stars, diff_aa_max_stars) among P/LP variants at the same codon."""
        parsed = _parse_codon(hgvsp)
        if parsed is None:
            return (0, 0)
        ref_aa, codon_num, alt_aa = parsed
        contig = self._resolve_contig(chrom)
        if contig is None:
            return (0, 0)
        p = int(pos)
        try:
            rows = list(self.tbx.fetch(contig, max(0, p - 4), p + 3))
        except ValueError:
            return (0, 0)
        same_aa_stars = 0
        diff_aa_stars = 0
        for line in rows:
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            r_pos = int(cols[1])
            r_ref = cols[3]
            r_alt = cols[4]
            if r_pos == p and r_ref == ref and r_alt == alt:
                continue
            info = _parse_vcf_info(cols[7])
            sig = _normalize_clnsig(info.get("CLNSIG", ""))
            if sig not in ("Pathogenic", "Likely_pathogenic"):
                continue
            stars = _clnrevstat_to_stars(info.get("CLNREVSTAT", ""))
            r_hgvsp = info.get("HGVSP") or info.get("CLNHGVS") or ""
            r_codon = _parse_codon(r_hgvsp)
            if r_codon is None:
                continue
            _, r_codon_num, r_alt_aa = r_codon
            if r_codon_num != codon_num:
                continue
            if r_alt_aa == alt_aa:
                same_aa_stars = max(same_aa_stars, stars)
            else:
                diff_aa_stars = max(diff_aa_stars, stars)
        return (same_aa_stars, diff_aa_stars)

    def same_codon(self, chrom: str, pos: int, ref: str, alt: str, hgvsp: str | None) -> tuple[bool, int] | None:
        """Return (same_aa: bool, stars: int) for a P/LP variant at the same codon, or None."""
        same_stars, diff_stars = self.codon_plp_stars(chrom, pos, ref, alt, hgvsp)
        if same_stars == 0 and diff_stars == 0:
            return None
        if same_stars >= diff_stars and same_stars > 0:
            return (True, same_stars)
        return (False, diff_stars)


def _parse_hotspots(raw: str) -> list[tuple[str, int, int]]:
    regions: list[tuple[str, int, int]] = []
    if not raw or not raw.strip():
        return regions
    for item in raw.replace(";", ",").split(","):
        item = item.strip()
        if not item or ":" not in item or "-" not in item:
            continue
        chrom_part, span = item.split(":", 1)
        start_s, end_s = span.split("-", 1)
        try:
            regions.append((norm_chrom(chrom_part), int(start_s), int(end_s)))
        except ValueError:
            continue
    return regions


def load_gene_table(path: Path | str) -> dict[str, dict]:
    """Load gene-level constraint, mechanism, ClinVar counts, hotspots, and BS1 thresholds."""
    p = Path(path)
    table: dict[str, dict] = {}
    with p.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue
            lof_raw = (row.get("lof_mechanism") or "").strip().lower()
            if lof_raw in ("yes", "true", "1"):
                lof_mech: bool | None = True
            elif lof_raw in ("no", "false", "0"):
                lof_mech = False
            else:
                lof_mech = None
            table[gene] = {
                "gene": gene,
                "loeuf": parse_float(row.get("loeuf")),
                "pli": parse_float(row.get("pli")),
                "moi": (row.get("moi") or "").strip() or None,
                "lof_mechanism": lof_mech,
                "clinvar_missense_plp": parse_int(row.get("clinvar_missense_plp")),
                "clinvar_missense_blb": parse_int(row.get("clinvar_missense_blb")),
                "clinvar_lof_plp": parse_int(row.get("clinvar_lof_plp")),
                "hotspot_regions": _parse_hotspots(row.get("hotspot_regions") or ""),
                "bs1_threshold": parse_float(row.get("bs1_threshold")),
            }
    return table


def read_slivar_info(vcf_path: Path | str) -> dict[tuple[str, int, str, str], dict]:
    """Extract gnomad_popmax_af and gnomad_nhomalt from a slivar-annotated VCF."""
    p = Path(vcf_path)
    out: dict[tuple[str, int, str, str], dict] = {}
    with pysam.VariantFile(str(p)) as vf:
        for rec in vf:
            chrom = norm_chrom(rec.chrom)
            pos = int(rec.pos)
            ref = str(rec.ref)
            popmax = parse_float(rec.info.get("gnomad_popmax_af"))
            nhomalt = parse_int(rec.info.get("gnomad_nhomalt"))
            for alt in (rec.alts or ()):
                out[(chrom, pos, ref, str(alt))] = {
                    "popmax_af": popmax,
                    "nhomalt": nhomalt,
                }
    return out
