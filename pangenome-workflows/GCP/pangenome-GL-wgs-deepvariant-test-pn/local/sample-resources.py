#!/usr/bin/env python3
"""Per-step container / process-tree resource sampler (wall time, peak RSS, cgroup max memory).

Supports both:
- Docker containers (cbontsr01): polls `docker ps` and reads `/sys/fs/cgroup/`
- Singularity / Apptainer under `cwltool --singularity` (IBEX): walks the process tree
  of each `singularity` / `Singularity runtime parent` step and sums RSS across all
  descendant processes in that step's subtree, so parallel CWL steps are measured
  independently.

Usage:
  python3 local/sample-resources.py --out /path/to/step-resources.tsv [--md /path/to/step-resources.md]
Send SIGTERM or SIGINT when the cwltool run finishes; it writes the final TSV and Markdown table on exit.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import signal
import subprocess
import sys
import time
from dataclasses import dataclass



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


def list_singularity_steps() -> list[tuple[str, str, str, int]]:
    """Return [(cid, image, command, subtree_rss_bytes), ...] for active Singularity/Apptainer steps."""
    try:
        out = subprocess.check_output(
            ["ps", "-eo", "pid=,ppid=,rss=,args="],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        return []

    procs: dict[int, tuple[int, int, str]] = {}
    children: dict[int, list[int]] = {}
    for line in out.splitlines():
        parts = line.strip().split(None, 3)
        if len(parts) < 4:
            continue
        try:
            pid = int(parts[0])
            ppid = int(parts[1])
            rss_kb = int(parts[2])
        except ValueError:
            continue
        args = parts[3]
        procs[pid] = (ppid, rss_kb * 1024, args)
        children.setdefault(ppid, []).append(pid)

    def is_sing_proc(args: str) -> bool:
        return (
            "Singularity runtime parent" in args
            or "Apptainer runtime parent" in args
            or args.startswith("singularity ")
            or args.startswith("apptainer ")
            or "/singularity " in args
            or "/apptainer " in args
        )

    # Restrict to descendants of our parent run.sh (if ppid > 1) so shared nodes
    # or overlapping runs never cross-contaminate step-resources.md
    parent_pid = os.getppid()
    allowed_pids: set[int] | None = None
    if parent_pid > 1 and parent_pid in procs:
        allowed_pids = set()
        q = [parent_pid]
        while q:
            cur = q.pop(0)
            if cur in allowed_pids:
                continue
            allowed_pids.add(cur)
            q.extend(children.get(cur, []))

    # Top-level singularity/apptainer roots (parent is not itself a singularity wrapper)
    roots: list[int] = []
    for pid, (ppid, _rss, args) in procs.items():
        if allowed_pids is not None and pid not in allowed_pids:
            continue
        if is_sing_proc(args):
            parent_args = procs.get(ppid, (0, 0, ""))[2]
            if not is_sing_proc(parent_args):
                roots.append(pid)

    res: list[tuple[str, str, str, int]] = []
    for root in roots:
        # BFS over subtree
        queue = [root]
        seen = set()
        total_rss = 0
        sif_name = ""
        first_payload_cmd = ""
        while queue:
            cur = queue.pop(0)
            if cur in seen or cur not in procs:
                continue
            seen.add(cur)
            _ppid, rss_b, c_args = procs[cur]
            total_rss += rss_b
            if not sif_name:
                m = re.search(r"([^\s/]+\.sif)\b", c_args)
                if m:
                    sif_name = m.group(1)
            if not is_sing_proc(c_args) and not first_payload_cmd:
                first_payload_cmd = c_args
            queue.extend(children.get(cur, []))
        root_args = procs[root][2]
        cmd = first_payload_cmd or root_args
        img = sif_name or (cmd.split()[0] if cmd else "singularity")
        res.append((f"sing-{root}", img, cmd, total_rss))
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

    # Ignore containers/steps that were already running before this workflow started
    baseline_cids = {cid for cid, _, _ in list_running_containers()} | {
        cid for cid, _, _, _ in list_singularity_steps()
    }
    records: dict[str, ContainerRecord] = {}
    stop = False

    def _handle_sig(_signum, _frame):
        nonlocal stop
        stop = True

    signal.signal(signal.SIGTERM, _handle_sig)
    signal.signal(signal.SIGINT, _handle_sig)

    while not stop:
        now = time.time()
        docker_rows = list_running_containers()
        if docker_rows:
            for cid, image, cmd in docker_rows:
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
        for cid, image, cmd, subtree_rss in list_singularity_steps():
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
                if rec.image == "singularity" and image != "singularity":
                    rec.image = image
                if ("singularity" in rec.command or "runtime parent" in rec.command) and (
                    "singularity" not in cmd and "runtime parent" not in cmd
                ):
                    rec.command = cmd
            rec.peak_rss_bytes = max(rec.peak_rss_bytes, subtree_rss)
            rec.peak_cgroup_bytes = max(rec.peak_cgroup_bytes, subtree_rss)
            rec.samples += 1
        write_outputs(records, args.out, args.md)
        time.sleep(args.interval)

    write_outputs(records, args.out, args.md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
