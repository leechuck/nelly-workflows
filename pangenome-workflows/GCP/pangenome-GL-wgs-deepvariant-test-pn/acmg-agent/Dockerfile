# nelly-acmg-agent:v9 -- ACMG evaluation agent runtime + WeasyPrint PDF renderer.
#
# Build context is the acmg-agent/ directory:
#   docker build -t nelly-acmg-agent:v9 acmg-agent/
# (MODEL_URL defaults to bartowski's Qwen2.5-3B-Instruct Q4_K_M GGUF, ~2 GB;
# build on ws, not the laptop. In production the step defaults to the lab
# unimatrix01 Qwen3.8-27B endpoint via ACMG_LLM_BACKEND=openai, keeping the
# local 3B GGUF as the offline fallback inside the image.)
#
# WHY 3B, NOT 7B (design §8.4 named 7B as the default): the real-model smoke
# test on ws with 8 vcpus showed Qwen2.5-7B-Instruct Q4_K_M needs ~130 s for
# the FIRST tool round (prompt reprocessing of the pack dominates) and never
# reaches a submit inside the 6 min/variant / 12-call budget -- every variant
# timed out. Qwen2.5-3B-Instruct Q4_K_M reaches a valid submit_assessment
# round-trip in ~120-300 s under the same budgets. The 7B remains available
# via --build-arg MODEL_URL on bigger nodes (design's Qwen3-8B caveat: same
# problem observed with 7B, not tested further). See wp4-REPORT.md.
#
# Gotchas applied (notes-cwl-gotchas.md):
# - llama-cpp-python must be built from source: prebuilt CPU wheels are
#   musl-linked and fail on Debian images (gotcha 18) ->
#   build-essential + cmake and CMAKE_ARGS=-DGGML_NATIVE=OFF so the build
#   is portable across ws and the cluster CPU families.
# - NO ENTRYPOINT (gotcha 4): Arvados prepends the entrypoint to the CWL
#   command and would double baseCommand. CMD ["python"] only.
# - The model is downloaded at build time into /model/model.gguf so run
#   time needs no network (llama backend is fully local).
FROM python:3.11-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential cmake curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ENV CMAKE_ARGS=-DGGML_NATIVE=OFF

# Model BEFORE the pip layer is not needed; download after pip install but
# keep ARG declared once at image top so both stages see it.
ARG MODEL_URL=https://huggingface.co/bartowski/Qwen2.5-3B-Instruct-GGUF/resolve/main/Qwen2.5-3B-Instruct-Q4_K_M.gguf

RUN pip install --no-cache-dir llama-cpp-python==0.3.2 pysam

RUN mkdir -p /model \
    && curl -fL --retry 3 --retry-delay 5 "${MODEL_URL}" -o /model/model.gguf \
    && test "$(stat -c %s /model/model.gguf)" -gt 100000000

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libpango-1.0-0 libpangoft2-1.0-0 libharfbuzz0b libharfbuzz-subset0 fonts-dejavu \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir weasyprint pypdf

COPY acmg_agent /app/acmg_agent

ENV PYTHONPATH=/app \
    ACMG_MODEL_PATH=/model/model.gguf \
    ACMG_LLM_BACKEND=llama \
    ACMG_IMAGE_TAG=nelly-acmg-agent:v9

# No ENTRYPOINT. The CWL tool drives: `python -m acmg_agent.cli assess|extract`.
CMD ["python"]
