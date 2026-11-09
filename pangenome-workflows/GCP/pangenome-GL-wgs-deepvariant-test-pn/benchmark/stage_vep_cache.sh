#!/bin/bash
# Stage the Ensembl VEP 116 GRCh38 cache for the modernised pipeline (WP5).
#
# Runs on the workstation `ws` (docker + disk space live there, not on the
# laptop). Downloads the VEP 116 cache tarball from the Ensembl FTP, verifies
# its size against the server's Content-Length, unpacks it into
#   /mnt/data1/nelly-sim/vep-116/cache/homo_sapiens/116_GRCh38/
# (the tarball itself unpacks to homo_sapiens/116_GRCh38/...), and PRINTS the
# `arv-put` command for the supervisor to run. This script never runs arv-put
# itself and never touches an Arvados cluster.
#
# Usage (on ws):  bash stage_vep_cache.sh
# Resumable:      re-running skips a complete download, resumes a partial one.
set -euo pipefail

CACHE_URL_BASE="https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache"
TARBALL="homo_sapiens_vep_116_GRCh38.tar.gz"
CACHE_URL="$CACHE_URL_BASE/$TARBALL"
STAGE_DIR="/mnt/data1/nelly-sim/vep-116"
LOG_DIR="/mnt/data1/nelly-sim/logs"
TARGET="$STAGE_DIR/cache"

mkdir -p "$STAGE_DIR" "$LOG_DIR"
cd "$STAGE_DIR"

# 1. Confirm the exact filename still exists on the FTP server (listing first,
#    per WP5 brief) and read its Content-Length for the size check.
#    The directory listing is ~3 MB (one row per species), so it occasionally
#    comes back truncated -- retry, and accept a successful HEAD on the
#    tarball URL as confirmation if the listing keeps coming up short.
listing=""
for attempt in 1 2 3; do
    listing=$(curl -sS "$CACHE_URL_BASE/" || true)
    if echo "$listing" | grep -q "$TARBALL"; then
        break
    fi
    echo "listing attempt $attempt: name not seen (listing may be truncated); retrying" >&2
    sleep 2
    listing=""
done
if ! echo "$listing" | grep -q "$TARBALL"; then
    if curl -sIL "$CACHE_URL" | tr -d '\r' | grep -q '^HTTP/.* 200'; then
        echo "listing check inconclusive, but HEAD on $CACHE_URL returned 200 -- continuing"
    else
        echo "ERROR: $TARBALL not found (listing + HEAD both failed)" >&2
        exit 1
    fi
fi
# Resolve redirects (ftp.ensembl.org answers with a 302 to a mirror) before
# reading Content-Length, and fall back to the listing size (26G) if HEAD
# yields nothing useful.
expected_size=$(curl -sIL "$CACHE_URL" | tr -d '\r' \
    | awk 'tolower($1)=="content-length:" { size=$2 } END { print size+0 }')
if [ "$expected_size" -lt 1000000000 ]; then
    # ~26 GB expected; anything under 1 GB means we failed to read a length.
    echo "WARNING: could not determine Content-Length ($expected_size); falling back to listing size" >&2
    listing_size=$(echo "$listing" | tr -d '\r' \
        | grep -A0 "$TARBALL" | head -1)
    case "$listing_size" in
        *[G]) gb=${listing_size%[G]*}; gb=${gb##* }; expected_size=$(( ${gb%.*} * 1024 * 1024 * 1024 )) ;;
    esac
fi

# 2. Download (resumable) unless a complete tarball is already present.
actual_size=$(stat -c %s "$TARBALL" 2>/dev/null || echo 0)
if [ "$actual_size" -ge "$expected_size" ]; then
    echo "$TARBALL already complete ($actual_size bytes); skipping download"
else
    if [ "$actual_size" -gt 0 ]; then
        echo "resuming partial download ($actual_size / $expected_size bytes)"
    fi
    curl -sS -o "$TARBALL" -C - "$CACHE_URL"
    actual_size=$(stat -c %s "$TARBALL")
fi

# 3. Verify size.
if [ "$actual_size" -lt "$expected_size" ]; then
    echo "ERROR: downloaded $actual_size bytes, expected $expected_size" >&2
    echo "(re-run this script to resume the download)" >&2
    exit 1
fi
echo "size verified: $actual_size bytes"

# 4. Unpack into $TARGET so the layout is $TARGET/homo_sapiens/116_GRCh38/.
#    VEP then runs with --offline --dir_cache/-dir "$TARGET" (uploaded to Keep
#    as the `vep_dir` Directory input) --cache_version 116.
if [ ! -d "$TARGET/homo_sapiens/116_GRCh38" ]; then
    echo "unpacking into $TARGET ..."
    mkdir -p "$TARGET"
    tar -xzf "$TARBALL" -C "$TARGET"
fi
#    (`|| true`: with pipefail, `ls | head` dies of SIGPIPE and would abort
#    the script before the arv-put line below is printed.)
ls -la "$TARGET/homo_sapiens/116_GRCh38" | head || true

# 5. Print (do NOT run) the arv-put command.
echo
echo "=== upload to Keep (run by supervisor on ws) ==="
echo "arv-put --name vep-116-GRCh38-cache /mnt/data1/nelly-sim/vep-116/cache"
