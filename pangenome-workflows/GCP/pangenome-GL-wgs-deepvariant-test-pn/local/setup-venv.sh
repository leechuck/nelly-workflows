#!/usr/bin/env bash
# Set up the local cwltool + Arvados Python environment on ontolinator (cbontsr01).
# Usage: bash local/setup-venv.sh
set -euo pipefail

if [ -z "${NELLY_LOCAL_ROOT:-}" ]; then
    if [ -d /ibex/scratch ]; then
        WORK_ROOT="/ibex/scratch/projects/c2014/nelly-local"
    else
        WORK_ROOT="/data/nelly-local"
    fi
else
    WORK_ROOT="$NELLY_LOCAL_ROOT"
fi
mkdir -p "$WORK_ROOT"/{inputs,images,runs,tmp,logs,sif_cache}
export TMPDIR="$WORK_ROOT/tmp"
export UV_CACHE_DIR="$WORK_ROOT/tmp/uv-cache"
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

echo "[setup-venv] Installing cwltool 3.x, arvados-python-client, nodejs-wheel-binaries, pyyaml, pytest, pysam..."
uv pip install --python "$VENV_DIR/bin/python" \
    "cwltool>=3.1" \
    "arvados-python-client>=2.7" \
    nodejs-wheel-binaries \
    pyyaml \
    pytest \
    pysam

ln -sf ../lib/python3.11/site-packages/nodejs_wheel/bin/node "$VENV_DIR/bin/node"

echo "[setup-venv] Done:"
"$VENV_DIR/bin/python" --version
"$VENV_DIR/bin/cwltool" --version
"$VENV_DIR/bin/arv-ls" --version
"$VENV_DIR/bin/node" --version
