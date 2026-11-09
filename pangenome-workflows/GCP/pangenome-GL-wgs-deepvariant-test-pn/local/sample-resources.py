#!/usr/bin/env python3
"""Per-step Docker container resource sampler (wall time, peak RSS, cgroup max memory).

Polls running Docker containers every `interval` seconds and reads cgroup v1/v2 memory
stats directly from `/sys/fs/cgroup/` (falling back to `docker stats --no-stream`), so
short-lived and hours-long CWL step containers alike are recorded without needing
`/usr/bin/time -v` inside the container images.

Usage:
  python3 local/sample-resources.py --out /path/to/step-resources.tsv [--md /path/to/step-resources.md]
Send SIGTERM or SIGINT when the cwltool run finishes; it writes the final TSV and Markdown table on exit.
"""
from __future__ import annotations

import argparse
import glob
import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, field


@dataclass
class ContainerRecord:
    cid: str
    image: str
    command: str
    first_seen: float
    last_seen: float
    peak_rss_bytes: int = 0
    peak_cgroup_bytes: int = 0
    samples: int = 0


def read_cgroup_memory(cid: str) -> tuple[int, int]:
    """Return (rss_bytes, cgroup_max_bytes) for a container ID from cgroup v1 or v2."""
    rss = 0
    cg_max = 0
    # cgroup v1 (Ubuntu 20.04 on cbontsr01)
    v1_dirs = [
        f"/sys/fs/cgroup/memory/docker/{cid}",
        f"/sys/fs/cgroup/memory/system.slice/docker-{cid}.scope",
    ]
    for d in v1_dirs:
        if os.path.isdir(d):
            stat_path = os.path.join(d, "memory.stat")
            if os.path.exists(stat_path):
                try:
                    with open(stat_path) as fh:
                        for line in fh:
                            parts = line.strip().split()
                            if len(parts) == 2 and parts[0] == "total_rss":
                                rss = max(rss, int(parts[1]))
                except OSError:
                    pass
            max_path = os.path.join(d, "memory.max_usage_in_bytes")
            if os.path.exists(max_path):
                try:
                    with open(max_path) as fh:
                        cg_max = max(cg_max, int(fh.read().strip()))
                except (OSError, ValueError):
                    pass
            return rss, cg_max

    # cgroup v2 fallback
    for p in glob.glob(f"/sys/fs/cgroup/**/docker-{cid}.scope", recursive=True):
        stat_path = os.path.join(p, "memory.stat")
        if os.path.exists(stat_path):
            try:
                with open(stat_path) as fh:
                    for line in fh:
                        parts = line.strip().split()
                        if len(parts) == 2 and parts[0] == "anon":
                            rss = max(rss, int(parts[1]))
            except OSError:
                pass
        peak_path = os.path.join(p, "memory.peak")
        if os.path.exists(peak_path):
            try:
                with open(peak_path) as fh:
                    cg_max = max(cg_max, int(fh.read().strip()))
            except (OSError, ValueError):
                pass
    return rss, cg_max


def list_running_containers() -> list[tuple[str, str, str]]:
    """Return [(full_cid, image, command), ...] for currently running containers."""
    try:
        out = subprocess.check_output(
            ["docker", "ps", "--no-trunc", "--format", "{{.ID}}\t{{.Image}}\t{{.Command}}"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        return []
    res = []
    for line in out.splitlines():
        parts = line.split("\t", 2)
        if len(parts) == 3:
            res.append((parts[0], parts[1], parts[2]))
    return res


def write_outputs(records: dict[str, ContainerRecord], tsv_path: str, md_path: str | None) -> None:
    rows = sorted(records.values(), key=lambda r: r.first_seen)
    with open(tsv_path, "w") as fh:
        fh.write("container_id\timage\tcommand\twall_seconds\tpeak_rss_mib\tpeak_cgroup_mib\tsamples\n")
        for r in rows:
            wall = max(0.0, r.last_seen - r.first_seen)
            rss_mib = r.peak_rss_bytes / (1024 * 1024)
            cg_mib = r.peak_cgroup_bytes / (1024 * 1024)
            fh.write(
                f"{r.cid[:12]}\t{r.image}\t{r.command}\t{wall:.1f}\t{rss_mib:.1f}\t{cg_mib:.1f}\t{r.samples}\n"
            )
    if md_path:
        with open(md_path, "w") as fh:
            fh.write("| Container | Image | Command | Wall (s) | Peak RSS (MiB) | Peak Cgroup Mem incl. cache (MiB) |\n")
            fh.write("|---|---|---|---:|---:|---:|\n")
            for r in rows:
                wall = max(0.0, r.last_seen - r.first_seen)
                rss_mib = r.peak_rss_bytes / (1024 * 1024)
                cg_mib = r.peak_cgroup_bytes / (1024 * 1024)
                cmd_short = r.command.replace("|", "\\|")[:60]
                fh.write(
                    f"| `{r.cid[:12]}` | `{r.image}` | `{cmd_short}` | {wall:.1f} | {rss_mib:.1f} | {cg_mib:.1f} |\n"
                )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Output TSV file")
    ap.add_argument("--md", default=None, help="Optional output Markdown table")
    ap.add_argument("--interval", type=float, default=1.0, help="Polling interval in seconds")
    args = ap.parse_args()

    # Ignore containers that were already running before this workflow started (long-lived services on ontolinator)
    baseline_cids = {cid for cid, _, _ in list_running_containers()}
    records: dict[str, ContainerRecord] = {}
    stop = False

    def _handle_sig(_signum, _frame):
        nonlocal stop
        stop = True

    signal.signal(signal.SIGTERM, _handle_sig)
    signal.signal(signal.SIGINT, _handle_sig)

    while not stop:
        now = time.time()
        for cid, image, cmd in list_running_containers():
            if cid in baseline_cids:
                continue
            rec = records.get(cid)
            if rec is None:
                rec = ContainerRecord(
                    cid=cid,
                    image=image,
                    command=cmd,
                    first_seen=now,
                    last_seen=now,
                )
                records[cid] = rec
            else:
                rec.last_seen = now
            rss, cg_max = read_cgroup_memory(cid)
            rec.peak_rss_bytes = max(rec.peak_rss_bytes, rss)
            rec.peak_cgroup_bytes = max(rec.peak_cgroup_bytes, cg_max)
            rec.samples += 1
        write_outputs(records, args.out, args.md)
        time.sleep(args.interval)

    write_outputs(records, args.out, args.md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
