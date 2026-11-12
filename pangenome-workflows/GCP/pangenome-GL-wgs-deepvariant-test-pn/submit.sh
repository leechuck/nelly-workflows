#!/usr/bin/env bash
# Submit a workflow to an Arvados cluster with a dev/prod intermediate-output policy.
#
#   ./submit.sh dev  <workflow.cwl> <inputs.json> [name]   # keep intermediate step
#                                                          # outputs for debugging (default)
#   ./submit.sh prod <workflow.cwl> <inputs.json> [name]   # trash intermediate step
#                                                          # outputs on SUCCESS (declared
#                                                          # workflow outputs are kept)
#
# Cluster: ARVADOS_API_HOST / ARVADOS_API_TOKEN from the environment; if unset, the
# KAUST cborg cluster and ~/.config/arvados/settings.conf are used.
#
# LLM key: committed input JSONs carry "llm_api_key": "FROM-ENV". This script replaces
# it with $ACMG_LLM_API_KEY (or the contents of ~/.config/borg-llm/api-key) in a
# mode-600 temporary copy, submits that copy and deletes it. The key never appears on
# a command line. main-vg.cwl declares llm_api_key a cwltool:Secrets input, so Arvados
# keeps it in secret_mounts rather than in the visible run inputs.
#
# Rationale for dev/prod: each step writes its own Keep collection. During development we
# want every intermediate; in production only the declared outputs, cleaned up only if
# the run succeeded (--trash-intermediate), so a failed run stays debuggable.
set -euo pipefail

MODE="${1:-dev}"; WF="${2:-}"; INPUTS="${3:-}"; NAME="${4:-}"
if [ -z "$WF" ] || [ -z "$INPUTS" ]; then
  echo "usage: $0 <dev|prod> <workflow.cwl> <inputs.json> [name]" >&2; exit 2
fi

case "$MODE" in
  dev)  TRASH=--no-trash-intermediate ;;
  prod) TRASH=--trash-intermediate ;;
  *) echo "mode must be 'dev' or 'prod', got '$MODE'" >&2; exit 2 ;;
esac
[ -z "$NAME" ] && NAME="$(basename "$WF" .cwl) ($MODE)"

if [ -z "${ARVADOS_API_HOST:-}" ]; then
  export ARVADOS_API_HOST=cborg.projectnelly.com
fi
if [ -z "${ARVADOS_API_TOKEN:-}" ]; then
  ARVADOS_API_TOKEN=$(python3 -c "import configparser,os; c=configparser.ConfigParser(); c.read_string('[d]\n'+open(os.path.expanduser('~/.config/arvados/settings.conf')).read()); print(c['d']['ARVADOS_API_TOKEN'])")
  export ARVADOS_API_TOKEN
fi

# Replace the FROM-ENV key placeholder in a private temporary copy of the inputs.
SUBMIT_INPUTS="$INPUTS"
TMP_INPUTS=""
cleanup() { [ -n "$TMP_INPUTS" ] && rm -f "$TMP_INPUTS"; return 0; }
trap cleanup EXIT
NEEDS_KEY=$(python3 - "$INPUTS" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
backend = str(d.get("llm_backend", "openai")).strip().lower()
print("yes" if d.get("llm_api_key") == "FROM-ENV" and backend == "openai" else
      "placeholder-unused" if d.get("llm_api_key") == "FROM-ENV" else "no")
PY
)
if [ "$NEEDS_KEY" != "no" ]; then
  KEY_FILE="$HOME/.config/borg-llm/api-key"
  if [ -z "${ACMG_LLM_API_KEY:-}" ] && [ -r "$KEY_FILE" ]; then
    ACMG_LLM_API_KEY=$(tr -d '[:space:]' < "$KEY_FILE")
  fi
  if [ -z "${ACMG_LLM_API_KEY:-}" ] && [ "$NEEDS_KEY" = "yes" ]; then
    echo "submit: $INPUTS needs an LLM key (llm_api_key is FROM-ENV, llm_backend is openai)." >&2
    echo "        Set ACMG_LLM_API_KEY or create $KEY_FILE, or set llm_backend to \"off\"." >&2
    exit 2
  fi
  umask 077
  TMP_INPUTS=$(mktemp "${TMPDIR:-$HOME/tmp}/submit-inputs.XXXXXX")
  mv "$TMP_INPUTS" "$TMP_INPUTS.json"; TMP_INPUTS="$TMP_INPUTS.json"
  ACMG_LLM_API_KEY="${ACMG_LLM_API_KEY:-}" python3 - "$INPUTS" "$TMP_INPUTS" <<'PY'
import json, os, sys
src, dst = sys.argv[1], sys.argv[2]
d = json.load(open(src))
base = os.path.dirname(os.path.abspath(src))
def absolutise(v):
    # Relative File/Directory paths in the original JSON resolve against its folder.
    if isinstance(v, dict):
        if v.get("class") in ("File", "Directory"):
            for k in ("path", "location"):
                p = v.get(k)
                if isinstance(p, str) and "://" not in p and not p.startswith(("keep:", "/")):
                    v[k] = os.path.join(base, p)
        for x in v.values():
            absolutise(x)
    elif isinstance(v, list):
        for x in v:
            absolutise(x)
absolutise(d)
d["llm_api_key"] = os.environ["ACMG_LLM_API_KEY"]
with open(dst, "w") as fh:
    json.dump(d, fh, indent=2)
PY
  SUBMIT_INPUTS="$TMP_INPUTS"
fi

echo "submit: host=$ARVADOS_API_HOST mode=$MODE ($TRASH) name='$NAME'" >&2
# Not exec: the temporary inputs copy must be removed after the runner has uploaded it.
set +e
arvados-cwl-runner --submit --no-wait "$TRASH" --name "$NAME" "$WF" "$SUBMIT_INPUTS"
RC=$?
set -e
exit $RC
