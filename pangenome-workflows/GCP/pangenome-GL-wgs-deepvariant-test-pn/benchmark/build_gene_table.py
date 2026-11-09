#!/usr/bin/env python3
"""Build the ACMG gene table (Supervisor: bake as `acmg_gene_table`, a Keep input).

    build_gene_table.py \
        --clinvar-summary variant_summary.txt.gz \
        --gnomad-constraint gnomad.v4.1.constraint_metrics.tsv \
        --out gene_table.tsv

Columns the engine reads (acmg_agent.inputs.load_gene_table):
    gene loeuf pli moi lof_mechanism clinvar_missense_plp
    clinvar_missense_blb clinvar_lof_plp hotspot_regions bs1_threshold

Per-gene ClinVar counts come from ClinVar `variant_summary.txt.gz` (one row
per allele): rows are grouped by the P/LP vs B/LB interpretation and the
Name field's HGVS (`c.` + `p.`) is classified as missense / LoF.

Consequence classes are approximations from the HGVS text only (variant_summary
carries no VEP-style consequence). Conservative rules:
    LoF   p.Ter / nonsense, p.? (incl. Met1?), any "fs", canonical splice site
          (c.N+/-1|2), counted also for start-loss;
    missense  single amino-acid substitution p.XaaNNNXaa (no Ter/fs/?/del/ins/dup);
    everything else is not counted.
MOI is left empty for every row (Exomiser supplies the disease MOI at run
time); hotspot_regions and bs1_threshold are left empty (no curated source
wired yet, engine treats them as unknown / uses defaults).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

OUT_COLUMNS = [
    "gene",
    "loeuf",
    "pli",
    "moi",
    "lof_mechanism",
    "clinvar_missense_plp",
    "clinvar_missense_blb",
    "clinvar_lof_plp",
    "hotspot_regions",
    "bs1_threshold",
]

LOEUF_DEFAULT_MAX = 0.35    # design-lite LOEUF threshold for LoF mechanism
MIN_CLEAN_COUNT = 10        # missense composition thresholds

_CODON_SUB_RE = re.compile(r"^p\.\(?([A-Za-z][a-z]{2})(\d+)\)?([A-Za-z][a-z]{2})?$")
_SPLICE_RE = re.compile(r"c\.\d+[+-][12](?!\d)")


def _col(row: dict, *names: str) -> str:
    """First present column value (gnomAD v4 renamed columns across releases)."""
    for n in names:
        v = row.get(n)
        if v is not None and str(v).strip() not in ("", "NA", "."):
            return str(v).strip()
    return ""


def classify_significance(raw: str) -> str | None:
    """'Pathogenic/Likely pathogenic, risk factor' -> 'PLP'; benign group -> 'BLB'."""
    segs = set()
    for part in raw.split("/"):
        part = part.split(",")[0].strip().lower()
        if part:
            segs.add(part)
    if "pathogenic" in segs or "likely_pathogenic" in segs or "likely pathogenic" in segs:
        return "PLP"
    if "benign" in segs or "likely_benign" in segs or "likely benign" in segs:
        return "BLB"
    return None


def classify_hgvs(name: str) -> str | None:
    """Approximate missense / LoF class from the variant_summary Name HGVS."""
    pm = re.search(r"p\.[^)\s]*", name)
    ptext = pm.group(0) if pm else ""
    cm = re.search(r"c\.[^)\s]*", name)
    ctext = cm.group(0) if cm else ""

    # canonical splice site (donor +1/+2, acceptor -1/-2)
    if _SPLICE_RE.search(ctext):
        return "LOF"
    # nonsense / frameshift / unknown protein consequence / start loss
    if ptext and ("Ter" in ptext or "fs" in ptext or "?" in ptext):
        return "LOF"
    # clean single-residue substitution -> missense
    if _CODON_SUB_RE.match(ptext):
        return "MISSENSE"
    return None


def read_gene_table_rows(summary_path: Path | str, assembly: str = "GRCh38") -> dict[str, dict]:
    """Aggregate per-gene missense/LoF P/LP and B/LB allele counts."""
    counts: dict[tuple[str, str, str], int] = {}
    opener = gzip.open if str(summary_path).endswith(".gz") else open
    with opener(summary_path, "rt", newline="") as f:  # type: ignore[operator]
        reader = csv.DictReader((line for line in f if not line.startswith("##")), delimiter="\t")
        header = reader.fieldnames or []
        header_map = {h.lstrip("#"): h for h in header}
        assembly_col = header_map.get("Assembly")
        for row in reader:
            sym = _col(row, "GeneSymbol")
            if not sym:
                continue
            if assembly and assembly_col and _col(row, "Assembly") != assembly:
                continue
            sig = classify_significance(_col(row, "ClinicalSignificance"))
            if sig is None:
                continue
            cls = classify_hgvs(_col(row, "Name"))
            if cls is None:
                continue
            key = (sym, sig, cls)
            counts.setdefault(key, 0)
            counts[key] += 1
    genes: dict[str, dict[str, int]] = {}
    for (sym, sig, cls), n in counts.items():
        gene = genes.setdefault(sym, {
            "clinvar_missense_plp": 0,
            "clinvar_missense_blb": 0,
            "clinvar_lof_plp": 0,
        })
        if sig == "PLP":
            if cls == "MISSENSE":
                gene["clinvar_missense_plp"] = n
            elif cls == "LOF":
                gene["clinvar_lof_plp"] = n
        elif sig == "BLB" and cls == "MISSENSE":
            gene["clinvar_missense_blb"] = n
    return genes


def read_constraint(path: Path | str) -> dict[str, dict[str, float | None]]:
    """gene -> {loeuf, pli} from gnomAD constraint metrics (v4.1 columns).

    gnomAD v4.1 lists one row per transcript; keep the MANE Select / canonical
    transcript's metrics for the gene, else the first row with values, and
    fill any NA fields from later rows of the same gene.
    """
    out: dict[str, dict[str, float | None]] = {}
    pref: dict[str, tuple[bool, bool]] = {}   # gene -> (mane, canonical) flags
    try:
        with open(path, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                gene = _col(row, "gene", "gene_symbol", "symbol")
                if not gene:
                    continue
                loeuf = _parse_opt_float(_col(row, "loeuf", "lof.oe_ci.upper", "lof.oe"))
                pli = _parse_opt_float(_col(row, "pli", "lof.pLI", "lof_hc_lc.pLI"))
                is_mane = _col(row, "mane_select").lower() == "true"
                is_canon = _col(row, "canonical").lower() == "true"
                cur = out.get(gene)
                if cur is None:
                    out[gene] = {"loeuf": loeuf, "pli": pli}
                    pref[gene] = (is_mane, is_canon)
                    continue
                if (is_mane or is_canon) and not pref[gene][0] and not pref[gene][1]:
                    # upgrade to a preferred transcript's values
                    if loeuf is not None:
                        cur["loeuf"] = loeuf
                    if pli is not None:
                        cur["pli"] = pli
                    pref[gene] = (is_mane, is_canon)
                else:
                    if cur["loeuf"] is None:
                        cur["loeuf"] = loeuf
                    if cur["pli"] is None:
                        cur["pli"] = pli
    except OSError as exc:
        raise RuntimeError(f"cannot open gnomAD constraint file {path}: {exc}") from exc
    return out


def _parse_opt_float(raw: str) -> float | None:
    try:
        return float(raw)
    except ValueError:
        return None


def lof_mechanism(loeuf: float | None, lof_plp: int | None, missense_plp: int | None) -> str:
    """yes / no / unknown rule from the WP4 brief."""
    lof = lof_plp or 0
    mis = missense_plp or 0
    if lof >= 10 or (loeuf is not None and loeuf < LOEUF_DEFAULT_MAX):
        return "yes"
    if lof == 0 and mis >= MIN_CLEAN_COUNT:
        return "no"
    return "unknown"


def build(out_path: Path | str, clinvar_summary: Path | str,
          gnomad_constraint: Path | str, assembly: str = "GRCh38") -> int:
    genes = read_gene_table_rows(clinvar_summary, assembly=assembly)
    constraint = read_constraint(gnomad_constraint)
    all_genes = sorted(set(genes) | set(constraint))
    written = 0
    try:
        with open(out_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(OUT_COLUMNS)
            for gene in all_genes:
                g = genes.get(gene, {})
                c = constraint.get(gene, {})
                loeuf = c.get("loeuf")
                pli = c.get("pli")
                mis_plp = g.get("clinvar_missense_plp")
                mis_blb = g.get("clinvar_missense_blb")
                lof_plp = g.get("clinvar_lof_plp")
                writer.writerow([
                    gene,
                    "" if loeuf is None else f"{loeuf:.6g}",
                    "" if pli is None else f"{pli:.6g}",
                    "",  # moi: Exomiser supplies the disease MOI at run time
                    lof_mechanism(loeuf, lof_plp, mis_plp),
                    "" if mis_plp is None else mis_plp,
                    "" if mis_blb is None else mis_blb,
                    "" if lof_plp is None else lof_plp,
                    "",  # hotspot_regions: no curated source wired
                    "",  # bs1_threshold: engine defaults (AD 1e-4 / AR 5e-3)
                ])
                written += 1
    except OSError as exc:
        raise RuntimeError(f"cannot write gene table {out_path}: {exc}") from exc
    return written


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clinvar-summary", required=True,
                   help="ClinVar variant_summary.txt[.gz]")
    p.add_argument("--gnomad-constraint", required=True,
                   help="gnomAD constraint metrics TSV (v4.1 style columns)")
    p.add_argument("--assembly", default="GRCh38",
                   help="variant_summary Assembly filter (default GRCh38)")
    p.add_argument("--out", required=True, help="output gene_table.tsv path")
    args = p.parse_args(argv)
    _ = Path(args.clinvar_summary).exists(), Path(args.gnomad_constraint).exists()
    n = build(args.out, args.clinvar_summary, args.gnomad_constraint, assembly=args.assembly)
    print(f"wrote {n} gene rows to {args.out}")
    sys.exit(0)


if __name__ == "__main__":
    main()
