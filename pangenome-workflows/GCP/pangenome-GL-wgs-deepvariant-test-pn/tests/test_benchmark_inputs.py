"""Guard that every GIAB benchmark input JSON supplies all of the workflow's REQUIRED inputs.

Extends the pattern from test_inputs_complete.py to benchmark-giab.cwl and
benchmark-giab-legacy.cwl against their corresponding JSON input files.
"""
import glob
import json
import os

import pytest
import yaml

WF_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCHMARK_DIR = os.path.join(WF_DIR, "benchmark")
BENCHMARK_CWL = os.path.join(BENCHMARK_DIR, "benchmark-giab.cwl")
LEGACY_CWL = os.path.join(BENCHMARK_DIR, "benchmark-giab-legacy.cwl")

BENCHMARK_JSONS = sorted(glob.glob(os.path.join(BENCHMARK_DIR, "benchmark-giab.HG*.json")))
LEGACY_JSONS = sorted(glob.glob(os.path.join(BENCHMARK_DIR, "benchmark-giab-legacy*.json")))


def _is_optional(type_val):
    """A CWL input type is optional if it is nullable (`T?`, or a list containing 'null')."""
    if isinstance(type_val, str):
        return type_val.endswith("?")
    if isinstance(type_val, list):
        return "null" in type_val
    return False


def required_inputs(cwl_path):
    with open(cwl_path) as fh:
        doc = yaml.safe_load(fh)
    required = set()
    for name, spec in (doc.get("inputs") or {}).items():
        if isinstance(spec, dict):
            if "default" in spec:  # has a default -> optional
                continue
            if _is_optional(spec.get("type")):
                continue
            required.add(name)
        else:  # bare-string shorthand, e.g. "File" / "File?"
            if not _is_optional(spec):
                required.add(name)
    return required


REQUIRED_BENCHMARK = required_inputs(BENCHMARK_CWL)
REQUIRED_LEGACY = required_inputs(LEGACY_CWL)


def test_benchmark_has_required_inputs():
    expected = {"graph", "hapl", "ref", "reads1", "reads2", "dv_mode", "truth_vcf", "confident_bed"}
    for name in expected:
        assert name in REQUIRED_BENCHMARK, f"{name} should be a required workflow input of benchmark-giab.cwl"


def test_legacy_has_required_inputs():
    expected = {"graph", "dist", "min", "ref", "reads1", "reads2", "truth_vcf", "confident_bed"}
    for name in expected:
        assert name in REQUIRED_LEGACY, f"{name} should be a required workflow input of benchmark-giab-legacy.cwl"


def test_benchmark_jsons_exist():
    assert len(BENCHMARK_JSONS) == 4, (
        f"Expected exactly 4 benchmark-giab JSON files (HG002/HG003 x standard/pangenome_aware), "
        f"found {len(BENCHMARK_JSONS)}: {BENCHMARK_JSONS}"
    )


def test_legacy_jsons_exist():
    assert len(LEGACY_JSONS) >= 1, "No benchmark-giab-legacy*.json files found"


@pytest.mark.parametrize("path", BENCHMARK_JSONS, ids=lambda p: os.path.basename(p))
def test_benchmark_json_provides_all_required_inputs(path):
    with open(path) as fh:
        provided = json.load(fh)
    missing = sorted(r for r in REQUIRED_BENCHMARK if provided.get(r) is None)
    assert not missing, (
        f"{os.path.basename(path)} is missing required workflow input(s): {missing}. "
        f"Every required input of benchmark-giab.cwl must be supplied or the run cannot launch."
    )


@pytest.mark.parametrize("path", LEGACY_JSONS, ids=lambda p: os.path.basename(p))
def test_legacy_json_provides_all_required_inputs(path):
    with open(path) as fh:
        provided = json.load(fh)
    missing = sorted(r for r in REQUIRED_LEGACY if provided.get(r) is None)
    assert not missing, (
        f"{os.path.basename(path)} is missing required workflow input(s): {missing}. "
        f"Every required input of benchmark-giab-legacy.cwl must be supplied."
    )
