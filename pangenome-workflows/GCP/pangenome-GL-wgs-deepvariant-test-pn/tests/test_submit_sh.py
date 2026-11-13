"""R12 Task 4: submit.sh replaces the FROM-ENV key placeholder privately and never puts the key on argv."""
from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
SUBMIT = HERE / "submit.sh"
KEY = "sk-unit-test-key-0001"

STUB = """#!/usr/bin/env bash
printf '%s\\n' "$@" > "$STUB_DIR/argv.txt"
last="${@: -1}"
cp "$last" "$STUB_DIR/inputs.json"
stat -c %a "$last" > "$STUB_DIR/mode.txt"
echo "$last" > "$STUB_DIR/path.txt"
exit 0
"""


def _setup(tmp_path, inputs):
    stub_dir = tmp_path / "stub"
    stub_dir.mkdir()
    runner = stub_dir / "arvados-cwl-runner"
    runner.write_text(STUB)
    runner.chmod(0o755)
    home = tmp_path / "home"
    (home / "tmp").mkdir(parents=True)
    inp = tmp_path / "inputs.json"
    inp.write_text(json.dumps(inputs))
    env = {
        "PATH": f"{stub_dir}:{os.environ['PATH']}",
        "HOME": str(home),
        "TMPDIR": str(home / "tmp"),
        "STUB_DIR": str(stub_dir),
        "ARVADOS_API_HOST": "arvados.example.org",
        "ARVADOS_API_TOKEN": "token-not-real",
    }
    return stub_dir, inp, env, home


def _run(env, inp):
    return subprocess.run(["bash", str(SUBMIT), "prod", "main-vg.cwl", str(inp)],
                          env=env, capture_output=True, text=True)


def test_key_injected_private_and_cleaned(tmp_path):
    stub, inp, env, home = _setup(tmp_path, {"llm_api_key": "FROM-ENV", "llm_backend": "openai",
                                              "reads1": {"class": "File", "path": "r1.fq.gz"}})
    env["ACMG_LLM_API_KEY"] = KEY
    r = _run(env, inp)
    assert r.returncode == 0, r.stderr
    sent = json.loads((stub / "inputs.json").read_text())
    assert sent["llm_api_key"] == KEY
    assert sent["reads1"]["path"] == str(tmp_path / "r1.fq.gz")  # relative path kept resolvable
    assert (stub / "mode.txt").read_text().strip() == "600"
    assert KEY not in (stub / "argv.txt").read_text()
    assert KEY not in r.stdout + r.stderr
    assert "--trash-intermediate" in (stub / "argv.txt").read_text()
    assert not Path((stub / "path.txt").read_text().strip()).exists()  # temp copy removed
    assert json.loads(inp.read_text())["llm_api_key"] == "FROM-ENV"  # original untouched
    assert "arvados.example.org" in r.stderr


def test_key_from_file(tmp_path):
    stub, inp, env, home = _setup(tmp_path, {"llm_api_key": "FROM-ENV"})
    (home / ".config" / "borg-llm").mkdir(parents=True)
    (home / ".config" / "borg-llm" / "api-key").write_text(KEY + "\n")
    r = _run(env, inp)
    assert r.returncode == 0, r.stderr
    assert json.loads((stub / "inputs.json").read_text())["llm_api_key"] == KEY


def test_missing_key_exits_2(tmp_path):
    stub, inp, env, home = _setup(tmp_path, {"llm_api_key": "FROM-ENV"})
    r = _run(env, inp)
    assert r.returncode == 2
    assert "needs an LLM key" in r.stderr
    assert not (stub / "argv.txt").exists()


def test_llm_off_needs_no_key(tmp_path):
    stub, inp, env, home = _setup(tmp_path, {"llm_api_key": "FROM-ENV", "llm_backend": "off"})
    r = _run(env, inp)
    assert r.returncode == 0, r.stderr


def test_no_placeholder_passes_file_unchanged(tmp_path):
    stub, inp, env, home = _setup(tmp_path, {"llm_backend": "off"})
    r = _run(env, inp)
    assert r.returncode == 0, r.stderr
    assert (stub / "path.txt").read_text().strip() == str(inp)


def test_endpoint_and_model_overrides_from_env(tmp_path):
    stub, inp, env, home = _setup(tmp_path, {"llm_api_key": "FROM-ENV", "llm_backend": "openai",
                                              "llm_endpoint": "http://lab.invalid:8000/v1",
                                              "llm_model": "qwen3.8-27b"})
    env.update(ACMG_LLM_API_KEY=KEY, ACMG_LLM_ENDPOINT="https://openrouter.ai/api/v1",
               ACMG_LLM_MODEL="qwen/qwen3.8-27b",
               ACMG_LLM_EXTRA_BODY='{"provider": {"data_collection": "deny"}}')
    r = _run(env, inp)
    assert r.returncode == 0, r.stderr
    sent = json.loads((stub / "inputs.json").read_text())
    assert sent["llm_endpoint"] == "https://openrouter.ai/api/v1"
    assert sent["llm_model"] == "qwen/qwen3.8-27b"
    assert json.loads(sent["llm_extra_body"]) == {"provider": {"data_collection": "deny"}}
    assert sent["llm_api_key"] == KEY
    assert KEY not in r.stdout + r.stderr
    assert json.loads(inp.read_text())["llm_model"] == "qwen3.8-27b"  # original untouched


def test_overrides_without_placeholder_keep_existing_key_field(tmp_path):
    stub, inp, env, home = _setup(tmp_path, {"llm_backend": "off"})
    env["ACMG_LLM_MODEL"] = "qwen/qwen3.8-27b"
    r = _run(env, inp)
    assert r.returncode == 0, r.stderr
    sent = json.loads((stub / "inputs.json").read_text())
    assert sent["llm_model"] == "qwen/qwen3.8-27b" and "llm_api_key" not in sent
