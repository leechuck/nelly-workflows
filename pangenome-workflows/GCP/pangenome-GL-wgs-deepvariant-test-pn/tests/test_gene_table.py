"""build_gene_table.py: columns and LoF-mechanism branches on a small fixture.

The fixture (tests/fixtures/gene_table/, 36 data lines covering HEXA, GENE_B
and GENE_C plus GRCh37/risk-factor/conflicting/no-symbol noise rows) exercises
the classifier paths offline; a run against the full real files is done on ws
(see wp4-REPORT.md).
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
FIX = HERE / "fixtures" / "gene_table"
SCRIPT = HERE.parent / "benchmark" / "build_gene_table.py"

EXPECTED_COLUMNS = [
    "gene", "loeuf", "pli", "moi", "lof_mechanism",
    "clinvar_missense_plp", "clinvar_missense_blb", "clinvar_lof_plp",
    "hotspot_regions", "bs1_threshold",
]


@pytest.fixture(scope="module")
def gene_table(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("gene_table") / "gene_table.tsv"
    res = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--clinvar-summary", str(FIX / "variant_summary.sub.txt"),
         "--gnomad-constraint", str(FIX / "gnomad.constraint.sub.tsv"),
         "--out", str(out)],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr
    return out


def _rows(path: Path) -> dict[str, dict]:
    with path.open() as f:
        return {r["gene"]: r for r in csv.DictReader(f, delimiter="\t")}


def test_columns_exact(gene_table):
    with gene_table.open() as f:
        header = f.readline().rstrip("\n").split("\t")
    assert header == EXPECTED_COLUMNS


def test_hexa_row(gene_table):
    rows = _rows(gene_table)
    hexa = rows["HEXA"]
    assert hexa["moi"] == ""                # left empty; Exomiser supplies MOI
    assert hexa["bs1_threshold"] == ""
    assert hexa["hotspot_regions"] == ""
    assert hexa["lof_mechanism"] == "yes"   # LOEUF 0.30 < 0.35
    assert hexa["clinvar_missense_plp"] == "3"
    assert hexa["clinvar_missense_blb"] == "2"
    assert hexa["clinvar_lof_plp"] == "3"   # nonsense + p.? + splice, GRCh38 only


def test_lof_mechanism_no_branch(gene_table):
    rows = _rows(gene_table)
    # GENE_B: loeuf 1.2, no LoF P/LP, 15 missense P/LP -> "no"
    assert rows["GENE_B"]["lof_mechanism"] == "no"
    assert rows["GENE_B"]["clinvar_missense_plp"] == "15"
    assert rows["GENE_B"]["clinvar_lof_plp"] == "0"


def test_lof_mechanism_unknown_branch(gene_table):
    rows = _rows(gene_table)
    # GENE_C: 1 splice LoF P/LP + 5 missense P/LP, loeuf 0.7 -> "unknown"
    assert rows["GENE_C"]["lof_mechanism"] == "unknown"
    assert rows["GENE_C"]["clinvar_lof_plp"] == "1"


def test_grch37_and_conflicting_rows_excluded(gene_table):
    rows = _rows(gene_table)
    # the GRCh37 HEXA missense row would have inflated the count to 4
    assert rows["HEXA"]["clinvar_missense_plp"] == "3"


def test_engine_loads_table(gene_table):
    # cross-check through the real engine reader (subprocess keeps the
    # pipeline-level test tree free of an acmg_agent package dependency)
    code = (
        "import sys;sys.path.insert(0,sys.argv[1]);"
        "from acmg_agent.inputs import load_gene_table;"
        "t=load_gene_table(sys.argv[2]);"
        "print(repr({k:(t[k]['loeuf'],t[k]['lof_mechanism'],t[k]['clinvar_lof_plp'],t[k]['bs1_threshold']) for k in t}))"
    )
    res = subprocess.run(
        [sys.executable, "-c", code, str(HERE.parent / "acmg-agent"), str(gene_table)],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr
    import ast

    data = ast.literal_eval(res.stdout)
    assert data["HEXA"][0] == pytest.approx(0.3)
    assert data["HEXA"][1] is True                       # "yes"
    assert data["HEXA"][2] == 3
    assert data["HEXA"][3] is None                       # empty -> engine default
    assert data["GENE_B"][1] is False                    # "no"
    assert data["GENE_C"][1] is None                     # "unknown"
