#!/usr/bin/env bash
# Set up the local cwltool + Arvados Python environment on ontolinator (cbontsr01).
# Usage: bash local/setup-venv.sh
set -euo pipefail

WORK_ROOT="${NELLY_LOCAL_ROOT:-/data/nelly-local}"
mkdir -p "$WORK_ROOT"/{inputs,images,runs,tmp,logs}
export TMPDIR="$WORK_ROOT/tmp"
export UV_LINK_MODE=copy

export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
if ! command -v uv >/dev/null 2>&1; then
    echo "[setup-venv] Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
fi

VENV_DIR="$WORK_ROOT/venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "[setup-venv] Creating Python 3.11 venv at $VENV_DIR..."
    uv venv --python 3.11 "$VENV_DIR"
fi

echo "[setup-venv] Installing cwltool 3.x, arvados-python-client, pyyaml, pytest, pysam..."
uv pip install --python "$VENV_DIR/bin/python" \
    "cwltool>=3.1" \
    "arvados-python-client>=2.7" \
    pyyaml \
    pytest \
    pysam

echo "[setup-venv] Done:"
"$VENV_DIR/bin/python" --version
"$VENV_DIR/bin/cwltool" --version
"$VENV_DIR/bin/arv-ls" --version
