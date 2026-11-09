#!/usr/bin/env bash
# Wrapper to run a CWL workflow locally on ontolinator (cbontsr01) with cwltool:
# - Exports TMPDIR=/data/nelly-local/tmp (never /tmp)
# - Injects "llm_api_key" from ~/.config/borg-llm/api-key into a runtime-only mode-600
#   copy of the inputs JSON when the template specifies "llm_api_key": "FROM-ENV"
#   (the key is never committed or logged)
# - Runs local/sample-resources.py alongside cwltool to record per-step container
#   wall time and peak RSS into <outdir>/step-resources.{tsv,md}
#
# Usage:
#   bash local/run.sh <outdir> <workflow.cwl> <inputs.json> [extra cwltool args...]
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
WORK_ROOT="${NELLY_LOCAL_ROOT:-/data/nelly-local}"
export TMPDIR="$WORK_ROOT/tmp"
mkdir -p "$TMPDIR" "$OUTDIR"

PYTHON="$WORK_ROOT/venv/bin/python"
CWLTOOL="$WORK_ROOT/venv/bin/cwltool"

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

"$PYTHON" - "$INPUTS_JSON" "$RUNTIME_JSON" <<'PY'
import json, os, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src) as fh:
    data = json.load(fh)
if data.get("llm_api_key") == "FROM-ENV":
    key_path = os.path.expanduser("~/.config/borg-llm/api-key")
    if os.path.exists(key_path):
        with open(key_path) as kf:
            data["llm_api_key"] = kf.read().strip()
    else:
        data["llm_api_key"] = os.environ.get("ACMG_LLM_API_KEY", "")
# Resolve relative file/directory paths in the input JSON against the input JSON's directory
with open(dst, "w") as out:
    json.dump(data, out, indent=2)
PY

"$PYTHON" "$SCRIPT_DIR/sample-resources.py" \
  --out "$OUTDIR/step-resources.tsv" \
  --md "$OUTDIR/step-resources.md" &
SAMPLER_PID=$!

INHIBIT=()
if command -v systemd-inhibit >/dev/null 2>&1; then
  INHIBIT=(systemd-inhibit --what=idle:sleep:shutdown --who="nelly-local" --why="Running nelly-workflows CWL pipeline")
fi

echo "[$(date -Is)] Starting cwltool: workflow=$WORKFLOW inputs=$INPUTS_JSON outdir=$OUTDIR"
"${INHIBIT[@]}" "$CWLTOOL" \
  --parallel \
  --timestamps \
  --tmpdir-prefix "$TMPDIR/" \
  --tmp-outdir-prefix "$TMPDIR/out-" \
  --outdir "$OUTDIR" \
  "$@" \
  "$WORKFLOW" \
  "$RUNTIME_JSON"
RC=$?
echo "[$(date -Is)] cwltool finished with exit code $RC"
exit $RC
