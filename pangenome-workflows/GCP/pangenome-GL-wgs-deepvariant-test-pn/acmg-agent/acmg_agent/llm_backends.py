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

import http.client
import json
import os
from pathlib import Path
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Protocol
from urllib.parse import urlparse


@dataclass
class ToolCall:
    """One function call requested by the model."""
    name: str
    arguments: dict
    id: str | None = None


@dataclass
class Message:
    """Role/content message; tool results carry ``role="tool"`` + ``name``.

    ``tool_calls`` is non-empty only on assistant messages that request
    function calls. ``finish_reason`` and ``usage`` capture completion
    metadata from the backend response when present.
    """
    role: str
    content: str | None = None
    tool_calls: list[ToolCall] = field(default_factory=list)
    name: str | None = None
    tool_call_id: str | None = None
    finish_reason: str | None = None
    usage: dict[str, Any] | None = None
    _trace_logged: bool = field(default=False, repr=False, compare=False)


def append_llm_trace(
    entry: dict[str, Any],
    trace_path: str | os.PathLike[str] | None = None,
) -> None:
    """Append one JSONL trace record to ``<outdir>/llm-trace.jsonl``.

    Never logs API keys or Authorization headers. Errors during trace
    writing are swallowed so diagnostics can never crash an assessment.
    """
    raw_path = trace_path or os.environ.get("ACMG_LLM_TRACE_PATH", "").strip()
    if not raw_path:
        return
    try:
        p = Path(raw_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, default=str) + "\n")
    except OSError:
        pass


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
                "id": tc.id or f"call_{idx}",
                "type": "function",
                "function": {
                    "name": tc.name,
                    "arguments": json.dumps(tc.arguments),
                },
            }
            for idx, tc in enumerate(msg.tool_calls)
        ]
    elif msg.role == "tool":
        out["content"] = msg.content or ""
        if msg.name:
            out["name"] = msg.name
        out["tool_call_id"] = msg.tool_call_id or msg.name or "call_0"
    return out


def _parse_tool_calls(raw_tool_calls: list[dict] | None) -> list[ToolCall]:
    calls: list[ToolCall] = []
    for tc in raw_tool_calls or []:
        fn = (tc.get("function") or {}) if isinstance(tc, dict) else {}
        name = (fn.get("name") or "").strip()
        if not name:
            continue
        tc_id = tc.get("id") if isinstance(tc, dict) else None
        if tc_id is not None and not isinstance(tc_id, str):
            tc_id = str(tc_id)
        args_raw = fn.get("arguments")
        try:
            args = json.loads(args_raw) if isinstance(args_raw, str) else (args_raw or {})
        except (json.JSONDecodeError, TypeError):
            # Invalid JSON: surface as an empty call; the agent's
            # contract enforcement will reject it and the model may retry.
            args = {"__invalid_json__": args_raw if args_raw is not None else ""}
        if not isinstance(args, dict):
            args = {"__invalid_args__": str(args)}
        calls.append(ToolCall(name=name, arguments=args, id=tc_id))
    return calls


def _completion_to_message(resp: dict) -> Message:
    choice = (resp.get("choices") or [{}])[0]
    m = choice.get("message") or {}
    # Explicitly ignore m.get("reasoning") and m.get("reasoning_content"):
    # chain-of-thought is never fed back as assistant content.
    content = m.get("content")
    if content is not None and not isinstance(content, str):
        content = str(content)
    finish_reason = choice.get("finish_reason")
    if finish_reason is not None and not isinstance(finish_reason, str):
        finish_reason = str(finish_reason)
    usage = resp.get("usage") if isinstance(resp.get("usage"), dict) else None
    return Message(
        role="assistant",
        content=content,
        tool_calls=_parse_tool_calls(m.get("tool_calls")),
        finish_reason=finish_reason,
        usage=usage,
    )


def _is_truncated_empty(resp: dict, msg: Message) -> bool:
    choice = (resp.get("choices") or [{}])[0]
    finish_reason = choice.get("finish_reason")
    has_content = bool(msg.content and msg.content.strip())
    has_tools = bool(msg.tool_calls)
    has_broken_tool_json = any(
        isinstance(tc.arguments, dict) and "__invalid_json__" in tc.arguments
        for tc in (msg.tool_calls or [])
    )
    return finish_reason == "length" and ((not has_content and not has_tools) or has_broken_tool_json)


def _build_trace_entry(
    *,
    backend_name: str,
    model_name: str | None,
    session_affinity: str | None,
    url: str | None,
    max_tokens: int | None,
    n_messages: int,
    n_tools: int,
    resp: dict | None,
    msg: Message,
    elapsed_seconds: float,
    trace_meta: dict[str, Any] | None = None,
    attempt: str = "initial",
) -> dict[str, Any]:
    choice = ((resp or {}).get("choices") or [{}])[0] if isinstance(resp, dict) else {}
    choice_msg = (choice.get("message") or {}) if isinstance(choice, dict) else {}
    reasoning_text = choice_msg.get("reasoning_content") or choice_msg.get("reasoning") or ""
    meta = trace_meta or {}
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "backend": backend_name,
        "model": model_name,
        "session_affinity": session_affinity,
        "gene": meta.get("gene"),
        "phase": meta.get("phase", "chat"),
        "attempt": attempt,
        "request": {
            "url": url,
            "max_tokens": max_tokens,
            "n_messages": n_messages,
            "n_tools": n_tools,
        },
        "response": {
            "id": (resp or {}).get("id") if isinstance(resp, dict) else None,
            "model": (resp or {}).get("model") if isinstance(resp, dict) else model_name,
            # Routers such as OpenRouter name the upstream provider that served the call.
            "provider": (resp or {}).get("provider") if isinstance(resp, dict) else None,
            "finish_reason": msg.finish_reason,
            "usage": msg.usage,
            "reasoning_chars": len(str(reasoning_text)) if reasoning_text else 0,
            "content": msg.content,
            "content_chars": len(msg.content or ""),
            "tool_calls": [
                {"id": tc.id, "name": tc.name, "arguments": tc.arguments}
                for tc in (msg.tool_calls or [])
            ],
        },
        "elapsed_seconds": round(elapsed_seconds, 3),
    }


def _max_tokens_from_env(default: int = 4096) -> int:
    raw = os.environ.get("ACMG_LLM_MAX_TOKENS", "").strip()
    if not raw:
        return default
    try:
        val = int(raw)
        return val if val > 0 else default
    except ValueError:
        return default


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
        # 16384: 12 tool rounds of pack-sized tool messages overflow 8192
        # (real-model smoke: 'Requested tokens (8193) exceed context window');
        # models are trained on 32k, CPU KV cache stays small.
        n_ctx: int = 16384,
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
        self.trace_path: str | None = None

    def chat(
        self,
        messages: list[Message],
        tools: list[dict],
        *,
        max_tokens: int | None = None,
        trace_meta: dict[str, Any] | None = None,
    ) -> Message:
        eff_max_tokens = max_tokens if (max_tokens is not None and max_tokens > 0) else self.max_tokens
        payload = [_message_to_openai(m) for m in messages]
        kwargs: dict[str, Any] = {
            "messages": payload,
            "max_tokens": eff_max_tokens,
        }
        if tools:
            kwargs["tools"] = tools
            kwargs["tool_choice"] = "auto"
        t0 = time.monotonic()
        resp = self._llm.create_chat_completion(**kwargs)
        msg = _completion_to_message(resp)
        msg._trace_logged = True
        append_llm_trace(
            _build_trace_entry(
                backend_name=self.backend_name,
                model_name=self.model_name,
                session_affinity=None,
                url=None,
                max_tokens=eff_max_tokens,
                n_messages=len(payload),
                n_tools=len(tools or []),
                resp=resp,
                msg=msg,
                elapsed_seconds=time.monotonic() - t0,
                trace_meta=trace_meta,
            ),
            trace_path=self.trace_path,
        )
        return msg


KEY_PLACEHOLDER = "FROM-ENV"


def _extra_body_from_env() -> dict[str, Any]:
    """Optional provider-specific request fields from ACMG_LLM_EXTRA_BODY (a JSON object).

    Example for OpenRouter: {"provider": {"data_collection": "deny"}}.
    Anything that is not a JSON object is a configuration error.
    """
    raw = os.environ.get("ACMG_LLM_EXTRA_BODY", "").strip()
    if not raw:
        return {}
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"ACMG_LLM_EXTRA_BODY is not valid JSON: {exc}") from None
    if not isinstance(value, dict):
        raise ValueError("ACMG_LLM_EXTRA_BODY must be a JSON object")
    return value


def _connection_class(scheme: str):
    """HTTP(S) connection class for the endpoint (one place, patchable in tests)."""
    return http.client.HTTPSConnection if scheme == "https" else http.client.HTTPConnection


class LLMUnavailable(RuntimeError):
    """The LLM cannot be used for this run (key missing/rejected, endpoint unreachable)."""


class OpenAICompatBackend:
    """OpenAI-compatible endpoint (https or localhost server), stdlib only."""

    def __init__(
        self,
        endpoint: str,
        model: str,
        *,
        api_key_env: str = "ACMG_LLM_API_KEY",
        timeout_seconds: int = 300,
        max_tokens: int | None = None,
        max_retries: int = 2,
        backoff_seconds: float = 1.0,
        session_affinity: str | None = None,
        sample: str | None = None,
    ):
        endpoint = endpoint.rstrip("/")
        self.url = endpoint if endpoint.endswith("/chat/completions") else f"{endpoint}/chat/completions"
        # Only real chat-completion endpoints are contacted; never file:/custom schemes.
        if urlparse(self.url).scheme not in ("http", "https"):
            raise ValueError(f"ACMG_LLM_ENDPOINT must be an http(s) URL, got {self.url!r}")
        self.model = model
        self.api_key = os.environ.get(api_key_env, "")
        self.timeout_seconds = timeout_seconds
        self.max_tokens = max_tokens if max_tokens is not None else _max_tokens_from_env(4096)
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds
        sample_id = (sample or os.environ.get("ACMG_SAMPLE") or "sample").strip() or "sample"
        self.session_affinity = session_affinity or f"{sample_id}-{uuid.uuid4()}"
        self.backend_name = "openai"
        self.model_name = model
        self.trace_path: str | None = None
        self.extra_body = _extra_body_from_env()

    def preflight(self, timeout_seconds: int = 30) -> None:
        """Fail fast with a plain reason before any review starts.

        Raises LLMUnavailable when the key is missing (or still the FROM-ENV
        placeholder), when the endpoint rejects the key, or when it cannot be
        reached. Never includes the key in the message.
        """
        if not self.api_key.strip() or self.api_key.strip() == KEY_PLACEHOLDER:
            raise LLMUnavailable(
                "LLM API key not supplied (the FROM-ENV placeholder was not replaced at submission)"
            )
        base = self.url[: -len("/chat/completions")] if self.url.endswith("/chat/completions") else self.url
        parsed = urlparse(base + "/models")
        conn_cls = _connection_class(parsed.scheme)
        conn = conn_cls(parsed.netloc, timeout=timeout_seconds)
        try:
            conn.request("GET", parsed.path, headers={"Authorization": f"Bearer {self.api_key}"})
            status = int(getattr(conn.getresponse(), "status", 200))
        except (OSError, http.client.HTTPException) as exc:
            reason = str(exc).replace(self.api_key, "***")[:120] or exc.__class__.__name__
            raise LLMUnavailable(f"LLM endpoint unreachable: {reason}") from None
        finally:
            conn.close()
        if status in (401, 403):
            raise LLMUnavailable(f"LLM endpoint rejected the API key (HTTP {status})")
        if status >= 500:
            raise LLMUnavailable(f"LLM endpoint unreachable: HTTP {status}")
        self._preflight_chat(timeout_seconds=max(timeout_seconds, 60))

    def _preflight_chat(self, timeout_seconds: int = 60) -> None:
        """One tiny chat request: proves the key, the model name and the account.

        Needed because some providers (OpenRouter) serve GET /models without
        authentication, so a wrong key or model name would otherwise only show
        up as one failed review per variant. A slow or busy endpoint (timeout,
        429, 5xx) is not an error here; the review calls have their own retries.
        """
        body = dict(self.extra_body)
        body.update({
            "model": self.model,
            "messages": [{"role": "user", "content": "Reply with OK."}],
            "max_tokens": 16,
        })
        payload = json.dumps(body).encode()
        parsed = urlparse(self.url)
        conn = _connection_class(parsed.scheme)(parsed.netloc, timeout=timeout_seconds)
        try:
            conn.request("POST", parsed.path or "/", body=payload, headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {self.api_key}",
                "Content-Length": str(len(payload)),
            })
            resp = conn.getresponse()
            status = int(getattr(resp, "status", 200))
            detail = resp.read().decode("utf-8", errors="replace") if status >= 400 else ""
        except (OSError, http.client.HTTPException):
            return
        finally:
            conn.close()
        detail = " ".join(detail.replace(self.api_key, "***").split())[:160]
        if status in (401, 403):
            raise LLMUnavailable(f"LLM endpoint rejected the API key (HTTP {status})")
        if status == 402:
            raise LLMUnavailable("LLM account has no credit left (HTTP 402)")
        if status in (400, 404, 422):
            raise LLMUnavailable(
                f"LLM endpoint refused model '{self.model}' (HTTP {status}): {detail}"
            )

    def _post_json(self, body: dict[str, Any]) -> dict[str, Any]:
        payload = json.dumps(body).encode()
        parsed = urlparse(self.url)
        conn_cls = _connection_class(parsed.scheme)
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
            "Content-Length": str(len(payload)),
            "X-Session-Affinity": self.session_affinity,
        }
        last_err: Exception | None = None
        for attempt in range(self.max_retries + 1):
            conn = conn_cls(parsed.netloc, timeout=self.timeout_seconds)
            try:
                conn.request("POST", parsed.path or "/", body=payload, headers=headers)
                http_resp = conn.getresponse()
                status = int(getattr(http_resp, "status", 200))
                resp_data = http_resp.read().decode("utf-8", errors="replace")
            except (OSError, http.client.HTTPException) as exc:
                last_err = exc
                if attempt < self.max_retries:
                    if self.backoff_seconds > 0:
                        time.sleep(self.backoff_seconds * (2 ** attempt))
                    continue
                raise RuntimeError(f"LLM endpoint connection error after {attempt + 1} attempts: {exc}") from exc
            finally:
                conn.close()

            if status == 429 or 500 <= status < 600:
                last_err = RuntimeError(f"LLM endpoint returned HTTP {status}")
                if attempt < self.max_retries:
                    if self.backoff_seconds > 0:
                        time.sleep(self.backoff_seconds * (2 ** attempt))
                    continue
                raise last_err
            if status >= 400:
                raise RuntimeError(f"LLM endpoint returned HTTP {status}: {resp_data[:200]}")
            try:
                return json.loads(resp_data)
            except json.JSONDecodeError as e:
                raise RuntimeError(f"LLM endpoint returned non-JSON payload (HTTP {str(e)[:60]})") from e
        raise RuntimeError(f"LLM endpoint request failed: {last_err}")

    def chat(
        self,
        messages: list[Message],
        tools: list[dict],
        *,
        max_tokens: int | None = None,
        trace_meta: dict[str, Any] | None = None,
    ) -> Message:
        # http.client over a scheme-validated http(s) URL only (no file:/custom schemes).
        eff_max_tokens = max_tokens if (max_tokens is not None and max_tokens > 0) else self.max_tokens
        payload_msgs = [_message_to_openai(m) for m in messages]
        # Provider-specific options first, so they can never override the fields below.
        base_body: dict[str, Any] = dict(self.extra_body)
        base_body.update({
            "model": self.model,
            "messages": payload_msgs,
            "max_tokens": eff_max_tokens,
        })
        if tools:
            base_body["tools"] = tools
            base_body["tool_choice"] = "auto"

        t0 = time.monotonic()
        resp = self._post_json(base_body)
        msg = _completion_to_message(resp)
        msg._trace_logged = True
        append_llm_trace(
            _build_trace_entry(
                backend_name=self.backend_name,
                model_name=self.model_name,
                session_affinity=self.session_affinity,
                url=self.url,
                max_tokens=eff_max_tokens,
                n_messages=len(payload_msgs),
                n_tools=len(tools or []),
                resp=resp,
                msg=msg,
                elapsed_seconds=time.monotonic() - t0,
                trace_meta=trace_meta,
                attempt="initial",
            ),
            trace_path=self.trace_path,
        )
        if _is_truncated_empty(resp, msg):
            retry_max_tokens = eff_max_tokens * 2
            retry_msgs = list(payload_msgs) + [
                {
                    "role": "user",
                    "content": (
                        "Your previous reply exhausted max_tokens during internal reasoning "
                        "before emitting a tool call or response. Keep internal reasoning "
                        "brief (under 500 tokens) and emit your tool call or response directly."
                    ),
                }
            ]
            retry_body = dict(base_body, messages=retry_msgs, max_tokens=retry_max_tokens)
            t1 = time.monotonic()
            resp2 = self._post_json(retry_body)
            msg2 = _completion_to_message(resp2)
            msg2._trace_logged = True
            append_llm_trace(
                _build_trace_entry(
                    backend_name=self.backend_name,
                    model_name=self.model_name,
                    session_affinity=self.session_affinity,
                    url=self.url,
                    max_tokens=retry_max_tokens,
                    n_messages=len(payload_msgs),
                    n_tools=len(tools or []),
                    resp=resp2,
                    msg=msg2,
                    elapsed_seconds=time.monotonic() - t1,
                    trace_meta=trace_meta,
                    attempt="length_retry",
                ),
                trace_path=self.trace_path,
            )
            if _is_truncated_empty(resp2, msg2):
                raise RuntimeError(
                    f"LLM response truncated (finish_reason=length) with no content/tool_calls "
                    f"after retry at max_tokens={retry_max_tokens}"
                )
            return msg2
        return msg


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
            n_ctx = int(os.environ.get("ACMG_N_CTX", "16384"))
        except ValueError:
            n_ctx = 16384
        return LlamaCppBackend(path, n_ctx=n_ctx)
    if mode == "openai":
        endpoint = os.environ.get("ACMG_LLM_ENDPOINT", "").strip()
        model = os.environ.get("ACMG_LLM_MODEL", "").strip()
        if not endpoint or not model:
            return None
        return OpenAICompatBackend(endpoint, model)
    raise ValueError(f"unknown ACMG_LLM_BACKEND={mode!r} (expected llama|openai|off)")
