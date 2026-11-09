"""ClinGen Evidence Repository (ERepo) VCF builder and ACMG criterion-level concordance evaluator."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re
import sys
import tempfile
import pysam

EVALUATED_CRITERIA = (
    "PVS1", "PS1", "PM1", "PM2", "PM4", "PM5", "PP2", "PP3",
    "BA1", "BS1", "BS2", "BP1", "BP3", "BP4", "BP7",
)

CLASSES = ("Pathogenic", "Likely pathogenic", "VUS", "Likely benign", "Benign")

GRCH38_CONTIG_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569,
}

_NC_SNV_RE = re.compile(r"NC_0000(\d{2})\.(\d+):g\.(\d+)([ACGT]+)>([ACGT]+)")
_CODE_PREFIX_RE = re.compile(r"^(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])\b", re.IGNORECASE)
_VAR_STR_RE = re.compile(r"^(?:chr)?([0-9XYM]+):(\d+)\s+([ACGT]+)>([ACGT]+)$", re.IGNORECASE)


def _nc_num_to_chrom(num: int) -> str | None:
    if 1 <= num <= 22:
        return f"chr{num}"
    if num == 23:
        return "chrX"
    if num == 24:
        return "chrY"
    return None


def parse_nc_hgvs_from_expressions(expr_str: str) -> tuple[str, int, str, str] | None:
    """Extract GRCh38 (chrom, pos, ref, alt) from ERepo 'HGVS Expressions' string.

    When multiple assembly versions (e.g. .12 GRCh38, .11 GRCh37, .10 NCBI36) are present,
    selects the highest RefSeq version number (GRCh38). Skips non-SNV/MNV expressions.
    """
    if not expr_str:
        return None
    best: tuple[int, str, int, str, str] | None = None
    for m in _NC_SNV_RE.finditer(expr_str):
        chr_num = int(m.group(1))
        ver = int(m.group(2))
        pos = int(m.group(3))
        ref = m.group(4).upper()
        alt = m.group(5).upper()
        chrom = _nc_num_to_chrom(chr_num)
        if chrom is None:
            continue
        if best is None or ver > best[0]:
            best = (ver, chrom, pos, ref, alt)
    if best is None:
        return None
    _, chrom, pos, ref, alt = best
    return (chrom, pos, ref, alt)


def normalize_criterion_code(raw: str) -> str | None:
    """Normalize ERepo evidence code like 'PM2_Supporting' or 'PM3_Very Strong' to 'PM2' / 'PM3'."""
    s = (raw or "").strip()
    if not s:
        return None
    token = re.split(r"[_\-\s]", s, maxsplit=1)[0].upper()
    m = _CODE_PREFIX_RE.match(token)
    return m.group(1).upper() if m else None


def normalize_assertion(raw: str) -> str:
    low = (raw or "").strip().lower()
    if low == "pathogenic":
        return "Pathogenic"
    if low in ("likely pathogenic", "likely_pathogenic"):
        return "Likely pathogenic"
    if "uncertain" in low or low == "vus":
        return "VUS"
    if low in ("likely benign", "likely_benign"):
        return "Likely benign"
    if low == "benign":
        return "Benign"
    return raw.strip() or "VUS"


def parse_erepo_codes(raw_col: str) -> set[str]:
    codes: set[str] = set()
    for part in (raw_col or "").split(","):
        norm = normalize_criterion_code(part)
        if norm:
            codes.add(norm)
    return codes


def build_erepo_vcf(erepo_tsv: Path | str, out_vcf: Path | str) -> int:
    """Build a GRCh38 VCF (optionally bgzipped + tabix-indexed) from ERepo HGVS Expressions."""
    in_path = Path(erepo_tsv)
    out_path = Path(out_vcf)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    records: dict[tuple[str, int, str, str], tuple[str, str]] = {}
    with in_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            parsed = parse_nc_hgvs_from_expressions(row.get("HGVS Expressions") or "")
            if parsed is None:
                continue
            vid = (row.get("ClinVar Variation Id") or row.get("Allele Registry Id") or ".").strip() or "."
            gene = (row.get("HGNC Gene Symbol") or "").strip()
            records[parsed] = (vid, gene)

    chrom_order = list(GRCH38_CONTIG_LENGTHS.keys())
    sorted_keys = sorted(
        records.keys(),
        key=lambda k: (chrom_order.index(k[0]) if k[0] in chrom_order else 99, k[1], k[2], k[3]),
    )

    header_lines = ["##fileformat=VCFv4.2"]
    for cname, clen in GRCH38_CONTIG_LENGTHS.items():
        header_lines.append(f"##contig=<ID={cname},length={clen}>")
    header_lines.append('##INFO=<ID=GENE,Number=1,Type=String,Description="HGNC Gene Symbol">')
    header_lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")

    is_gz = str(out_path).endswith(".gz")
    raw_vcf_path = out_path.with_suffix("") if is_gz else out_path

    with raw_vcf_path.open("w") as out_f:
        for hl in header_lines:
            out_f.write(hl + "\n")
        for chrom, pos, ref, alt in sorted_keys:
            vid, gene = records[(chrom, pos, ref, alt)]
            info = f"GENE={gene}" if gene else "."
            out_f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/1\n")

    if is_gz:
        pysam.tabix_compress(str(raw_vcf_path), str(out_path), force=True)
        pysam.tabix_index(str(out_path), preset="vcf", force=True)
        if raw_vcf_path.exists():
            raw_vcf_path.unlink()

    return len(sorted_keys)


def _parse_variant_key(var_str: str) -> tuple[str, int, str, str] | None:
    m = _VAR_STR_RE.match((var_str or "").strip())
    if not m:
        return None
    chrom = m.group(1).upper()
    chrom = f"chr{chrom}" if chrom != "MT" else "chrM"
    return (chrom, int(m.group(2)), m.group(3).upper(), m.group(4).upper())


def compare_assessment_to_erepo(
    assessment_json: Path | str,
    erepo_tsv: Path | str,
    out_md: Path | str,
) -> dict:
    """Compare acmg-assessment.json against ERepo Met / Not Met columns and Assertions."""
    assessment = json.loads(Path(assessment_json).read_text())
    by_var: dict[tuple[str, int, str, str], dict] = {}
    for cand in assessment.get("candidates") or []:
        vkey = _parse_variant_key(cand.get("variant") or "")
        if vkey is not None:
            by_var[vkey] = cand

    per_crit: dict[str, dict[str, float | int]] = {
        code: {"tp": 0, "fp": 0, "fn": 0, "tn": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0}
        for code in EVALUATED_CRITERIA
    }
    class_matrix: Counter[tuple[str, str]] = Counter()
    matched = 0

    with Path(erepo_tsv).open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            vkey = parse_nc_hgvs_from_expressions(row.get("HGVS Expressions") or "")
            if vkey is None or vkey not in by_var:
                continue
            matched += 1
            cand = by_var[vkey]
            truth_met = parse_erepo_codes(row.get("Applied Evidence Codes (Met)") or "")
            truth_not_met = parse_erepo_codes(row.get("Applied Evidence Codes (Not Met)") or "")
            pred_by_code = {c["code"]: c["status"] for c in (cand.get("criteria") or [])}

            for code in EVALUATED_CRITERIA:
                pred_is_met = pred_by_code.get(code) == "met"
                if code in truth_met:
                    if pred_is_met:
                        per_crit[code]["tp"] += 1
                    else:
                        per_crit[code]["fn"] += 1
                else:
                    if pred_is_met:
                        per_crit[code]["fp"] += 1
                    elif code in truth_not_met:
                        per_crit[code]["tn"] += 1

            truth_cls = normalize_assertion(row.get("Assertion") or "")
            pred_cls = normalize_assertion(cand.get("classification") or "")
            class_matrix[(truth_cls, pred_cls)] += 1

    for code, d in per_crit.items():
        tp, fp, fn = int(d["tp"]), int(d["fp"]), int(d["fn"])
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0
        d["precision"] = round(prec, 4)
        d["recall"] = round(rec, 4)
        d["f1"] = round(f1, 4)

    lines = [
        "# ClinGen Evidence Repository (ERepo) Concordance Report",
        "",
        f"- **Matched variants evaluated**: {matched}",
        "",
        "## 1. Per-Criterion Concordance",
        "",
        "| Criterion | TP | FP | FN | TN | Precision | Recall | F1 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for code in EVALUATED_CRITERIA:
        d = per_crit[code]
        lines.append(
            f"| {code} | {d['tp']} | {d['fp']} | {d['fn']} | {d['tn']} | "
            f"{d['precision']:.3f} | {d['recall']:.3f} | {d['f1']:.3f} |"
        )

    lines.extend([
        "",
        "## 2. Final Classification Confusion Matrix (ERepo Truth Rows vs Engine Columns)",
        "",
        "| ERepo \\ Engine | " + " | ".join(CLASSES) + " |",
        "|---|" + "|".join(["---"] * len(CLASSES)) + "|",
    ])
    for t_cls in CLASSES:
        row_cells = [str(class_matrix[(t_cls, p_cls)]) for p_cls in CLASSES]
        lines.append(f"| **{t_cls}** | " + " | ".join(row_cells) + " |")

    out_path = Path(out_md)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")

    return {
        "matched_variants": matched,
        "per_criterion": per_crit,
        "class_matrix": class_matrix,
    }


def main(argv: list[str] | None = None) -> int:
    args_list = list(argv if argv is not None else sys.argv[1:])
    if args_list and args_list[0] in ("build-vcf", "vcf", "compare", "evaluate"):
        args_list = args_list[1:]

    p = argparse.ArgumentParser(description="ClinGen ERepo VCF builder and ACMG concordance evaluator")
    p.add_argument("--erepo", required=True, help="Path to ClinGen ERepo TSV")
    p.add_argument("--out-vcf", default=None, help="Output GRCh38 VCF (.vcf or .vcf.gz)")
    p.add_argument("--assessment", default=None, help="Path to acmg-assessment.json from assess --from-vcf")
    p.add_argument("--out", default=None, help="Output Markdown concordance report")
    ns = p.parse_args(args_list)

    if ns.out_vcf:
        build_erepo_vcf(ns.erepo, ns.out_vcf)
    if ns.assessment and ns.out:
        compare_assessment_to_erepo(ns.assessment, ns.erepo, ns.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
