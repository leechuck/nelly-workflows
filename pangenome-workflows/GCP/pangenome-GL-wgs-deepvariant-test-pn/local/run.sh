#!/usr/bin/env bash
# Wrapper to run a CWL workflow locally on ontolinator (cbontsr01) or IBEX (Slurm + Singularity) with cwltool:
# - Auto-detects HOST=cbontsr01|ibex (override with HOST=... or NELLY_LOCAL_ROOT=...)
# - Exports TMPDIR=$WORK_ROOT/tmp (never /tmp)
# - Rewrites /data/nelly-local paths in the input JSON to $WORK_ROOT when running on IBEX
# - Injects "llm_api_key" from ~/.config/borg-llm/api-key into a runtime-only mode-600
#   copy of the inputs JSON when the template specifies "llm_api_key": "FROM-ENV"
#   (the key is never committed or logged)
# - Runs local/sample-resources.py alongside cwltool to record per-step container / process-tree
#   wall time and peak RSS into <outdir>/step-resources.{tsv,md}
#
# Usage:
#   [HOST=cbontsr01|ibex] bash local/run.sh <outdir> <workflow.cwl> <inputs.json> [extra cwltool args...]
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <outdir> <workflow.cwl> <inputs.json> [extra cwltool args...]" >&2
  exit 2
fi

OUTDIR="$1"
WORKFLOW="$2"
INPUTS_JSON="$3"
shift 3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "${HOST:-}" ]; then
  if [ -d /ibex/scratch ] || [ -n "${SLURM_JOB_ID:-}" ]; then
    HOST="ibex"
  else
    HOST="cbontsr01"
  fi
fi

if [ "$HOST" = "ibex" ]; then
  WORK_ROOT="${NELLY_LOCAL_ROOT:-/ibex/scratch/projects/c2014/nelly-local}"
else
  WORK_ROOT="${NELLY_LOCAL_ROOT:-/data/nelly-local}"
fi

export TMPDIR="$WORK_ROOT/tmp"
mkdir -p "$TMPDIR" "$OUTDIR"

PYTHON="$WORK_ROOT/venv/bin/python"
CWLTOOL="$WORK_ROOT/venv/bin/cwltool"
export PATH="$WORK_ROOT/venv/bin:$PATH"

ENGINE_ARGS=()
if [ "$HOST" = "ibex" ]; then
  if ! command -v singularity >/dev/null 2>&1; then
    if [ -f /etc/profile.d/modules.sh ]; then
      # shellcheck disable=SC1091
      source /etc/profile.d/modules.sh
    fi
    module load singularity
  fi
  mkdir -p "$WORK_ROOT/sif_cache/.cache"
  export SINGULARITY_CACHEDIR="$WORK_ROOT/sif_cache/.cache"
  export CWL_SINGULARITY_CACHE="$WORK_ROOT/sif_cache"
  export SINGULARITY_TMPDIR="$TMPDIR"
  export SINGULARITY_BINDPATH="/ibex/scratch/projects/c2014,/encrypted/e3008,$TMPDIR"
  ENGINE_ARGS=(--singularity)
fi

RUNTIME_JSON="$(mktemp "$TMPDIR/cwl-inputs.XXXXXX.json")"
chmod 600 "$RUNTIME_JSON"
SAMPLER_PID=""

cleanup() {
  rm -f "$RUNTIME_JSON"
  if [ -n "$SAMPLER_PID" ] && kill -0 "$SAMPLER_PID" 2>/dev/null; then
    kill -TERM "$SAMPLER_PID" 2>/dev/null || true
    wait "$SAMPLER_PID" 2>/dev/null || true
  fi
}
trap cleanup EXIT

"$PYTHON" - "$INPUTS_JSON" "$RUNTIME_JSON" "$WORK_ROOT" <<'PY'
import json, os, socket, sys
src, dst, work_root = sys.argv[1], sys.argv[2], sys.argv[3]
with open(src) as fh:
    raw = fh.read()
if work_root != "/data/nelly-local":
    raw = raw.replace("/data/nelly-local", work_root)
data = json.loads(raw)
if data.get("llm_api_key") == "FROM-ENV":
    key_path = os.path.expanduser("~/.config/borg-llm/api-key")
    if os.path.exists(key_path):
        with open(key_path) as kf:
            data["llm_api_key"] = kf.read().strip()
    else:
        data["llm_api_key"] = os.environ.get("ACMG_LLM_API_KEY", "")
ep = data.get("llm_endpoint", "")
if "unimatrix01.kaust.edu.sa" in ep:
    has_ns = False
    try:
        with open("/etc/resolv.conf") as rf:
            has_ns = any(line.strip().startswith("nameserver") for line in rf)
        if has_ns:
            socket.gethostbyname("unimatrix01.kaust.edu.sa")
    except OSError:
        has_ns = False
    if not has_ns:
        data["llm_endpoint"] = ep.replace("unimatrix01.kaust.edu.sa", "10.72.186.139")
with open(dst, "w") as out:
    json.dump(data, out, indent=2)
PY

"$PYTHON" "$SCRIPT_DIR/sample-resources.py" \
  --out "$OUTDIR/step-resources.tsv" \
  --md "$OUTDIR/step-resources.md" &
SAMPLER_PID=$!

echo "[$(date -Is)] Starting cwltool (HOST=$HOST): workflow=$WORKFLOW inputs=$INPUTS_JSON outdir=$OUTDIR"
"$PYTHON" -c '
import sys, cwltool.job, cwltool.main
_orig_init = cwltool.job.JobBase.__init__
def _init_with_net(self, *args, **kwargs):
    _orig_init(self, *args, **kwargs)
    self.networkaccess = True
cwltool.job.JobBase.__init__ = _init_with_net
sys.exit(cwltool.main.run())
' \
  --parallel \
  --timestamps \
  "${ENGINE_ARGS[@]}" \
  --tmpdir-prefix "$TMPDIR/" \
  --tmp-outdir-prefix "$TMPDIR/out-" \
  --outdir "$OUTDIR" \
  "$@" \
  "$WORKFLOW" \
  "$RUNTIME_JSON"
RC=$?
echo "[$(date -Is)] cwltool finished with exit code $RC"
exit $RC
