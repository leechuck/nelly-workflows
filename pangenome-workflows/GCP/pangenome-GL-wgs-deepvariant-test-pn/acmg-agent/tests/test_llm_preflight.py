"""R12 Task 2: LLM failures (missing/rejected key, unreachable endpoint) are visible, fast, and never fail the step."""
from __future__ import annotations

import json
import socket
from pathlib import Path

import pytest

from acmg_agent import llm_backends
from acmg_agent.cli import main
from acmg_agent.llm_backends import LLMUnavailable, OpenAICompatBackend
from tests.test_pdf_report import _extract_pdf_text

FIX = Path(__file__).parent / "fixtures" / "acmg"
SECRET = "sk-test-secret-value-123"


class _Resp:
    def __init__(self, status):
        self.status = status

    def read(self):
        return b"{}"


def _fake_conn(status=None, exc=None, seen=None):
    class Conn:
        def __init__(self, netloc, timeout=None):
            pass

        def request(self, method, path, body=None, headers=None):
            if seen is not None:
                seen.append((method, path))
            if exc is not None:
                raise exc

        def getresponse(self):
            return _Resp(status)

        def close(self):
            pass

    return Conn


def _run(tmp_path, monkeypatch, key, conn):
    monkeypatch.setenv("ACMG_LLM_BACKEND", "openai")
    monkeypatch.setenv("ACMG_LLM_ENDPOINT", "http://llm.invalid:8000/v1")
    monkeypatch.setenv("ACMG_LLM_MODEL", "qwen3.8-27b")
    monkeypatch.setenv("ACMG_LLM_API_KEY", key)
    monkeypatch.setattr(llm_backends, "_connection_class", lambda scheme: conn)
    out = {k: tmp_path / v for k, v in
           (("md", "acmg-report.md"), ("json", "acmg-assessment.json"), ("pdf", "acmg-report.pdf"))}
    rc = main([
        "assess",
        "--genes-tsv", str(FIX / "TS2/exomiser.genes.tsv"),
        "--variants-tsv", str(FIX / "TS2/exomiser.variants.tsv"),
        "--vep-vcf", str(FIX / "vep.vcf"),
        "--clinvar-vcf", str(FIX / "clinvar.vcf.gz"),
        "--gene-table", str(FIX / "gene_table.tsv"),
        "--sample", "TS2",
        "--output-md", str(out["md"]),
        "--output-json", str(out["json"]),
        "--output-pdf", str(out["pdf"]),
    ])
    assert rc == 0
    data = json.loads(out["json"].read_text())
    _, pdf_text = _extract_pdf_text(out["pdf"])
    return data, out["md"].read_text(), pdf_text


@pytest.mark.parametrize("key", ["", "FROM-ENV"])
def test_missing_key_fails_visibly_without_http(tmp_path, monkeypatch, capsys, key):
    seen = []
    data, md, pdf = _run(tmp_path, monkeypatch, key, _fake_conn(status=200, seen=seen))
    assert seen == []  # no HTTP call at all
    assert data["llm"]["status"] == "failed"
    assert "FROM-ENV placeholder was not replaced" in data["llm"]["error"]
    assert "LLM API key not supplied" in md
    assert "LLM review not run" in pdf and "API key not supplied" in pdf
    assert "LLM API key not supplied" in capsys.readouterr().err
    # deterministic result still present
    assert data["leading_candidate"]["gene"] == "HEXA"
    assert data["leading_candidate"]["classification"] == "Pathogenic"


@pytest.mark.parametrize("status", [401, 403])
def test_rejected_key(tmp_path, monkeypatch, capsys, status):
    data, md, pdf = _run(tmp_path, monkeypatch, SECRET, _fake_conn(status=status))
    assert data["llm"]["status"] == "failed"
    assert data["llm"]["error"] == f"LLM endpoint rejected the API key (HTTP {status})"
    assert "rejected the API key" in pdf
    blob = json.dumps(data) + md + pdf + capsys.readouterr().err
    assert SECRET not in blob


def test_unreachable_endpoint_fails_fast(tmp_path, monkeypatch):
    seen = []
    data, md, pdf = _run(tmp_path, monkeypatch, SECRET,
                         _fake_conn(exc=socket.timeout("timed out"), seen=seen))
    assert len(seen) == 1  # one preflight attempt, no per-candidate retries
    assert data["llm"]["status"] == "failed"
    assert data["llm"]["error"].startswith("LLM endpoint unreachable:")
    assert "unreachable" in md and "unreachable" in pdf
    assert all(c.get("agent_status") in ("not_reviewed", None) for c in data["candidates"])


def test_preflight_ok_passes(monkeypatch):
    monkeypatch.setenv("ACMG_LLM_API_KEY", SECRET)
    monkeypatch.setattr(llm_backends, "_connection_class", lambda scheme: _fake_conn(status=200))
    OpenAICompatBackend("http://llm.invalid:8000/v1", "m").preflight()


def test_preflight_message_never_contains_key(monkeypatch):
    monkeypatch.setenv("ACMG_LLM_API_KEY", SECRET)
    monkeypatch.setattr(llm_backends, "_connection_class", lambda scheme: _fake_conn(exc=OSError(f"boom {SECRET}")))
    with pytest.raises(LLMUnavailable) as ei:
        OpenAICompatBackend("http://llm.invalid:8000/v1", "m").preflight()
    assert SECRET not in str(ei.value)


# --- R13: chat preflight (providers such as OpenRouter serve GET /models without a key) ---

def _routed_conn(chat_status=200, chat_body=b"{}", chat_exc=None, seen=None):
    class Conn:
        def __init__(self, netloc, timeout=None):
            self._method = None

        def request(self, method, path, body=None, headers=None):
            self._method = method
            if seen is not None:
                seen.append((method, path, json.loads(body) if body else None))
            if method == "POST" and chat_exc is not None:
                raise chat_exc

        def getresponse(self):
            resp = _Resp(200 if self._method == "GET" else chat_status)
            resp.read = lambda: chat_body
            return resp

        def close(self):
            pass

    return Conn


def _backend(monkeypatch, conn, model="qwen/qwen3.8-27b"):
    monkeypatch.setenv("ACMG_LLM_API_KEY", SECRET)
    monkeypatch.setattr(llm_backends, "_connection_class", lambda scheme: conn)
    return OpenAICompatBackend("https://openrouter.invalid/api/v1", model)


@pytest.mark.parametrize("status", [401, 403])
def test_chat_preflight_catches_key_rejected_after_public_models(monkeypatch, status):
    with pytest.raises(LLMUnavailable, match=f"rejected the API key \\(HTTP {status}\\)"):
        _backend(monkeypatch, _routed_conn(chat_status=status)).preflight()


def test_chat_preflight_names_unknown_model_and_masks_key(monkeypatch):
    body = f'{{"error": "{SECRET} is not a valid model ID"}}'.encode()
    with pytest.raises(LLMUnavailable) as ei:
        _backend(monkeypatch, _routed_conn(chat_status=404, chat_body=body), model="qwen3.8-27b").preflight()
    msg = str(ei.value)
    assert "refused model 'qwen3.8-27b' (HTTP 404)" in msg and "not a valid model ID" in msg
    assert SECRET not in msg


def test_chat_preflight_no_credit(monkeypatch):
    with pytest.raises(LLMUnavailable, match="no credit"):
        _backend(monkeypatch, _routed_conn(chat_status=402)).preflight()


@pytest.mark.parametrize("kw", [{"chat_exc": socket.timeout("timed out")}, {"chat_status": 429},
                                {"chat_status": 503}])
def test_chat_preflight_tolerates_a_busy_endpoint(monkeypatch, kw):
    _backend(monkeypatch, _routed_conn(**kw)).preflight()


def test_extra_body_sent_but_cannot_override_model(monkeypatch):
    seen = []
    monkeypatch.setenv("ACMG_LLM_EXTRA_BODY",
                       '{"provider": {"data_collection": "deny"}, "model": "other", "max_tokens": 1}')
    backend = _backend(monkeypatch, _routed_conn(seen=seen))
    backend.preflight()
    post = [b for m, _, b in seen if m == "POST"][0]
    assert post["provider"] == {"data_collection": "deny"}
    assert post["model"] == "qwen/qwen3.8-27b" and post["max_tokens"] == 16


@pytest.mark.parametrize("raw", ["not json", "[1, 2]"])
def test_bad_extra_body_is_a_clear_error(monkeypatch, raw):
    monkeypatch.setenv("ACMG_LLM_EXTRA_BODY", raw)
    with pytest.raises(ValueError, match="ACMG_LLM_EXTRA_BODY"):
        _backend(monkeypatch, _routed_conn())


def test_trace_records_upstream_provider():
    from acmg_agent.llm_backends import Message, _build_trace_entry
    entry = _build_trace_entry(
        backend_name="openai", model_name="qwen/qwen3.8-27b", session_affinity="s", url="u",
        max_tokens=16, n_messages=1, n_tools=0,
        resp={"id": "x", "model": "qwen/qwen3.8-27b", "provider": "SomeHost",
              "choices": [{"message": {"content": "OK"}, "finish_reason": "stop"}]},
        msg=Message(role="assistant", content="OK"), elapsed_seconds=0.1, trace_meta=None, attempt="initial",
    )
    assert entry["response"]["provider"] == "SomeHost"
