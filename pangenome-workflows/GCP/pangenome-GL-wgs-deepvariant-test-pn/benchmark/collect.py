#!/usr/bin/env bash
""":"
exec uv run python "$0" "$@"
"""
# Above header allows running directly or via python
import argparse
import csv
import json
import os
import sys
from typing import Any, Dict, List, Optional


YANG_REFERENCE_ROWS = [
    {
        "sample": "HG003 (Yang et al.)",
        "graph": "GRCh38",
        "caller": "BWA-MEM2 + DV 1.6",
        "snp_recall": "-",
        "snp_precision": "-",
        "snp_f1": "0.9956",
        "indel_recall": "-",
        "indel_precision": "-",
        "indel_f1": "-",
        "snp_errors": "29,052",
        "indel_errors": "7,064",
        "difficult_snp_recall": "0.964",
    },
    {
        "sample": "HG003 (Yang et al.)",
        "graph": "HPRC v1.0",
        "caller": "vg 1.54 + DV 1.6",
        "snp_recall": "-",
        "snp_precision": "-",
        "snp_f1": "0.9973",
        "indel_recall": "-",
        "indel_precision": "-",
        "indel_f1": "-",
        "snp_errors": "18,099",
        "indel_errors": "6,609",
        "difficult_snp_recall": "0.984",
    },
    {
        "sample": "HG003 (Yang et al.)",
        "graph": "Saudi (KSA)",
        "caller": "vg 1.54 + DV 1.6",
        "snp_recall": "-",
        "snp_precision": "-",
        "snp_f1": "0.9971",
        "indel_recall": "-",
        "indel_precision": "-",
        "indel_f1": "-",
        "snp_errors": "18,996",
        "indel_errors": "6,582",
        "difficult_snp_recall": "0.983",
    },
]


def parse_happy_summary_csv(csv_path: str) -> Dict[str, Dict[str, float]]:
    """Parse hap.py summary.csv rows for Type in {SNP, INDEL} with Filter == PASS."""
    metrics: Dict[str, Dict[str, float]] = {}
    with open(csv_path, mode="r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            var_type = row.get("Type", "").strip()
            filt = row.get("Filter", "").strip()
            if var_type in ("SNP", "INDEL") and filt == "PASS":
                try:
                    recall = float(row["METRIC.Recall"])
                    precision = float(row["METRIC.Precision"])
                    f1 = float(row["METRIC.F1_Score"])
                except (KeyError, ValueError) as exc:
                    raise ValueError(f"Error parsing metrics in {csv_path}: {exc}") from exc
                metrics[var_type] = {
                    "recall": recall,
                    "precision": precision,
                    "f1": f1,
                }
    return metrics


def normalize_caller(caller: str) -> str:
    """Normalize caller labels, ensuring legacy runs are clearly identified."""
    c_lower = caller.strip().lower()
    if c_lower in ("legacy", "vg1.54", "legacy_vg1.54"):
        return "legacy (vg 1.54 + DV 1.6)"
    return caller.strip()


def build_markdown(runs_data: List[Dict[str, Any]], include_yang_refs: bool = True) -> str:
    lines: List[str] = [
        "# GIAB Benchmark Evaluation Results",
        "",
        "## Performance Summary (Confident Regions, Filter = PASS)",
        "",
        "| sample | graph | caller | SNP recall | SNP precision | SNP F1 | INDEL recall | INDEL precision | INDEL F1 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for run in runs_data:
        sample = run.get("sample", "-")
        graph = run.get("graph", "-")
        caller = normalize_caller(run.get("caller", "-"))
        csv_path = run.get("summary_csv")

        if csv_path and os.path.isfile(csv_path):
            m = parse_happy_summary_csv(csv_path)
            snp = m.get("SNP", {})
            indel = m.get("INDEL", {})
            snp_recall = f"{snp['recall']:.4f}" if "recall" in snp else "-"
            snp_precision = f"{snp['precision']:.4f}" if "precision" in snp else "-"
            snp_f1 = f"{snp['f1']:.4f}" if "f1" in snp else "-"
            indel_recall = f"{indel['recall']:.4f}" if "recall" in indel else "-"
            indel_precision = f"{indel['precision']:.4f}" if "precision" in indel else "-"
            indel_f1 = f"{indel['f1']:.4f}" if "f1" in indel else "-"
        else:
            snp_recall = snp_precision = snp_f1 = "-"
            indel_recall = indel_precision = indel_f1 = "-"

        lines.append(
            f"| {sample} | {graph} | {caller} | {snp_recall} | {snp_precision} | {snp_f1} | "
            f"{indel_recall} | {indel_precision} | {indel_f1} |"
        )

    if include_yang_refs:
        lines.extend([
            "",
            "## Published Reference Baselines (Yang et al. — vg 1.54 + DV 1.6, HG003 vs GIAB v4.2.1)",
            "",
            "| Reference Study | Graph | Caller / Aligner | SNP F1 | SNP Errors | INDEL Errors | Difficult Regions SNP Recall |",
            "|---|---|---|---|---|---|---|",
        ])
        for y in YANG_REFERENCE_ROWS:
            lines.append(
                f"| {y['sample']} | {y['graph']} | {y['caller']} | {y['snp_f1']} | {y['snp_errors']} | "
                f"{y['indel_errors']} | {y['difficult_snp_recall']} |"
            )

    lines.append("")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Collect hap.py summary statistics into a markdown report.")
    parser.add_argument("--runs", required=True, help="Path to runs.json containing run specifications.")
    parser.add_argument("--out", default="results-2026-11.md", help="Output markdown filename.")
    args = parser.parse_args(argv)

    with open(args.runs, "r", encoding="utf-8") as fh:
        runs = json.load(fh)

    if not isinstance(runs, list):
        raise ValueError("runs.json must contain a list of objects.")

    md_content = build_markdown(runs)

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(md_content)

    print(f"Wrote benchmark results to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
