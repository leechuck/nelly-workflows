"""Tests for benchmark/collect.py hap.py metrics collection and reporting."""
import json
import os
import sys

import pytest

# Ensure pipeline dir is in sys.path so benchmark/ can be imported
WF_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if WF_DIR not in sys.path:
    sys.path.insert(0, WF_DIR)

from benchmark.collect import build_markdown, main, normalize_caller, parse_happy_summary_csv

FIXTURE_SUMMARY_CSV = os.path.join(WF_DIR, "tests", "fixtures", "happy", "summary.csv")


def test_parse_happy_summary_csv():
    metrics = parse_happy_summary_csv(FIXTURE_SUMMARY_CSV)

    assert "SNP" in metrics, "Expected SNP metrics in parsed summary"
    assert "INDEL" in metrics, "Expected INDEL metrics in parsed summary"

    snp = metrics["SNP"]
    assert snp["recall"] == pytest.approx(0.997567, rel=1e-5)
    assert snp["precision"] == pytest.approx(0.998123, rel=1e-5)
    assert snp["f1"] == pytest.approx(0.997845, rel=1e-5)

    indel = metrics["INDEL"]
    assert indel["recall"] == pytest.approx(0.991234, rel=1e-5)
    assert indel["precision"] == pytest.approx(0.993456, rel=1e-5)
    assert indel["f1"] == pytest.approx(0.992344, rel=1e-5)


def test_normalize_caller():
    assert normalize_caller("legacy") == "legacy (vg 1.54 + DV 1.6)"
    assert normalize_caller("LEGACY") == "legacy (vg 1.54 + DV 1.6)"
    assert normalize_caller("standard") == "standard"
    assert normalize_caller("pangenome_aware") == "pangenome_aware"


def test_build_markdown_formats_four_decimals():
    runs = [
        {
            "sample": "HG003",
            "graph": "KSA",
            "caller": "standard",
            "summary_csv": FIXTURE_SUMMARY_CSV,
        },
        {
            "sample": "HG003",
            "graph": "KSA",
            "caller": "legacy",
            "summary_csv": FIXTURE_SUMMARY_CSV,
        },
    ]

    md = build_markdown(runs)

    # Check 4 decimal formatting for metrics
    assert "0.9976" in md  # SNP recall rounded
    assert "0.9981" in md  # SNP precision rounded
    assert "0.9978" in md  # SNP F1 rounded
    assert "0.9912" in md  # INDEL recall rounded
    assert "0.9935" in md  # INDEL precision rounded
    assert "0.9923" in md  # INDEL F1 rounded

    # Check legacy labeling
    assert "legacy (vg 1.54 + DV 1.6)" in md

    # Check Yang reference rows
    assert "HG003 (Yang et al.)" in md
    assert "0.9956" in md
    assert "0.9973" in md
    assert "0.9971" in md
    assert "29,052" in md
    assert "18,099" in md
    assert "18,996" in md


def test_collect_cli_end_to_end(tmp_path):
    runs_file = tmp_path / "runs.json"
    out_file = tmp_path / "results.md"

    runs_data = [
        {
            "sample": "HG002",
            "graph": "KSA",
            "caller": "pangenome_aware",
            "summary_csv": FIXTURE_SUMMARY_CSV,
        }
    ]
    runs_file.write_text(json.dumps(runs_data), encoding="utf-8")

    rc = main(["--runs", str(runs_file), "--out", str(out_file)])
    assert rc == 0
    assert out_file.exists()

    content = out_file.read_text(encoding="utf-8")
    assert "HG002" in content
    assert "pangenome_aware" in content
    assert "0.9976" in content
