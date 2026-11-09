"""Pluggable LLM backends for the ACMG evaluation agent.

The design (spec §4.3 step 4) requires a pluggable backend with a minimal
chat interface:

    class Backend(Protocol):
        def chat(messages: list[Message], tools: list[dict]) -> Message

Three concrete implementations (or the absence of one -> deterministic
fall-back) are selected from the environment:

    ACMG_LLM_BACKEND = llama | openai | off
    ACMG_MODEL_PATH  = path to the local GGUF (llama backend)
    ACMG_LLM_ENDPOINT / ACMG_LLM_MODEL (openai backend, for hosted models)

The heavy ``llama_cpp`` import happens lazily inside ``LlamaCppBackend`` so
that the unit-test suite (FakeBackend only) never needs the package.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Protocol
from urllib.parse import urlparse


@dataclass
class ToolCall:
    """One function call requested by the model."""
    name: str
    arguments: dict


@dataclass
class Message:
    """Role/content message; tool results carry ``role="tool"`` + ``name``.

    ``tool_calls`` is non-empty only on assistant messages that request
    function calls.
    """
    role: str
    content: str | None = None
    tool_calls: list[ToolCall] = field(default_factory=list)
    name: str | None = None


class Backend(Protocol):
    """Chat protocol every backend must satisfy."""

    def chat(self, messages: list[Message], tools: list[dict]) -> Message:
        ...


def _message_to_openai(msg: Message) -> dict:
    out: dict[str, Any] = {"role": msg.role}
    if msg.content is not None:
        out["content"] = msg.content
    if msg.role == "assistant" and msg.tool_calls:
        out["content"] = msg.content or None
        out["tool_calls"] = [
            {
                "type": "function",
                "function": {
                    "name": tc.name,
                    "arguments": json.dumps(tc.arguments),
                },
            }
            for tc in msg.tool_calls
        ]
    elif msg.role == "tool":
        out["content"] = msg.content or ""
        if msg.name:
            out["name"] = msg.name
    return out


def _parse_tool_calls(raw_tool_calls: list[dict] | None) -> list[ToolCall]:
    calls: list[ToolCall] = []
    for tc in raw_tool_calls or []:
        fn = (tc.get("function") or {}) if isinstance(tc, dict) else {}
        name = (fn.get("name") or "").strip()
        if not name:
            continue
        args_raw = fn.get("arguments")
        try:
            args = json.loads(args_raw) if isinstance(args_raw, str) else (args_raw or {})
        except (json.JSONDecodeError, TypeError):
            # Invalid JSON: surface as an empty call; the agent's
            # contract enforcement will reject it and the model may retry.
            args = {"__invalid_json__": args_raw if args_raw is not None else ""}
        if not isinstance(args, dict):
            args = {"__invalid_args__": str(args)}
        calls.append(ToolCall(name=name, arguments=args))
    return calls


def _completion_to_message(resp: dict) -> Message:
    choice = (resp.get("choices") or [{}])[0]
    m = choice.get("message") or {}
    content = m.get("content")
    if content is not None and not isinstance(content, str):
        content = str(content)
    return Message(
        role="assistant",
        content=content,
        tool_calls=_parse_tool_calls(m.get("tool_calls")),
    )


class LlamaCppBackend:
    """Local GGUF backend via llama-cpp-python function calling.

    Built around ``create_chat_completion(messages, tools=..., tool_choice='auto')``
    with a function-calling chat format (default ``chatml-function-calling``;
    llama-cpp-python 0.3.2 supports it). ``n_ctx`` >= 8192 and ``n_threads``
    driven by ``ACMG_THREADS`` (the CWL layer passes ``$(runtime.cores)``).
    """

    def __init__(
        self,
        model_path: str,
        *,
        n_ctx: int = 8192,
        n_threads: int | None = None,
        chat_format: str = "chatml-function-calling",
        max_tokens: int = 2048,
    ):
        # Deferred import: the heavy native module only exists inside the
        # container image; unit tests use FakeBackend and never reach this path.
        import importlib
        llama_cpp = importlib.import_module("llama_cpp")

        threads = n_threads or _threads_from_env()
        self.model_path = model_path
        self.max_tokens = max_tokens
        self._llm = llama_cpp.Llama(
            model_path=model_path,
            n_ctx=n_ctx,
            n_threads=threads,
            n_gpu_layers=0,           # CPU-only sizing (spec: fits 8 cores / 16 GiB)
            chat_format=chat_format,
            verbose=False,
        )
        self.backend_name = "llama"
        self.model_name = os.path.basename(model_path)

    def chat(self, messages: list[Message], tools: list[dict]) -> Message:
        payload = [_message_to_openai(m) for m in messages]
        kwargs: dict[str, Any] = {
            "messages": payload,
            "max_tokens": self.max_tokens,
        }
        if tools:
            kwargs["tools"] = tools
            kwargs["tool_choice"] = "auto"
        resp = self._llm.create_chat_completion(**kwargs)
        return _completion_to_message(resp)


class OpenAICompatBackend:
    """OpenAI-compatible endpoint (https or localhost server), stdlib only."""

    def __init__(
        self,
        endpoint: str,
        model: str,
        *,
        api_key_env: str = "ACMG_LLM_API_KEY",
        timeout_seconds: int = 300,
        max_tokens: int = 2048,
    ):
        endpoint = endpoint.rstrip("/")
        self.url = endpoint if endpoint.endswith("/chat/completions") else f"{endpoint}/chat/completions"
        # Only real chat-completion endpoints are contacted; never file:/custom schemes.
        if urlparse(self.url).scheme not in ("http", "https"):
            raise ValueError(f"ACMG_LLM_ENDPOINT must be an http(s) URL, got {self.url!r}")
        self.model = model
        self.api_key = os.environ.get(api_key_env, "")
        self.timeout_seconds = timeout_seconds
        self.max_tokens = max_tokens
        self.backend_name = "openai"
        self.model_name = model

    def chat(self, messages: list[Message], tools: list[dict]) -> Message:
        # http.client over a scheme-validated http(s) URL only (no file:/custom schemes).
        import http.client

        body: dict[str, Any] = {
            "model": self.model,
            "messages": [_message_to_openai(m) for m in messages],
            "max_tokens": self.max_tokens,
        }
        if tools:
            body["tools"] = tools
            body["tool_choice"] = "auto"
        payload = json.dumps(body).encode()
        parsed = urlparse(self.url)
        conn_cls = http.client.HTTPSConnection if parsed.scheme == "https" else http.client.HTTPConnection
        conn = conn_cls(parsed.netloc, timeout=self.timeout_seconds)
        try:
            conn.request("POST", parsed.path or "/", body=payload, headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {self.api_key}",
                "Content-Length": str(len(payload)),
            })
            resp_data = conn.getresponse().read().decode()
        finally:
            conn.close()
        try:
            resp = json.loads(resp_data)
        except json.JSONDecodeError as e:
            raise RuntimeError(f"LLM endpoint returned non-JSON payload (HTTP {str(e)[:60]})") from e
        return _completion_to_message(resp)


def _threads_from_env() -> int:
    try:
        return int(os.environ.get("ACMG_THREADS", os.cpu_count() or 1))
    except ValueError:
        return os.cpu_count() or 1


def select_backend_from_env() -> Backend | None:
    """Build a backend from ACMG_* environment variables; None = LLM off.

    ``ACMG_LLM_BACKEND`` selects ``llama`` (default), ``openai`` or ``off``.
    The llama path needs ``ACMG_MODEL_PATH`` (baked into the container image
    at /model/model.gguf); without it we return None instead of crashing so
    the deterministic path still runs.
    """
    mode = (os.environ.get("ACMG_LLM_BACKEND") or "llama").strip().lower()
    if mode in ("off", "none", "0", "disabled", ""):
        return None
    if mode == "llama":
        path = os.environ.get("ACMG_MODEL_PATH", "").strip()
        if not path:
            return None
        try:
            n_ctx = int(os.environ.get("ACMG_N_CTX", "8192"))
        except ValueError:
            n_ctx = 8192
        return LlamaCppBackend(path, n_ctx=n_ctx)
    if mode == "openai":
        endpoint = os.environ.get("ACMG_LLM_ENDPOINT", "").strip()
        model = os.environ.get("ACMG_LLM_MODEL", "").strip()
        if not endpoint or not model:
            return None
        return OpenAICompatBackend(endpoint, model)
    raise ValueError(f"unknown ACMG_LLM_BACKEND={mode!r} (expected llama|openai|off)")
