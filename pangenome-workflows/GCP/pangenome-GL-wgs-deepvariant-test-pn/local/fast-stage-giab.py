#!/usr/bin/env python3
"""16-worker parallel Arvados Keep block streamer + Singularity legacy SIF pre-puller.

Coordinates safely with local/sbatch-benchmark-giab.sh via <dest>.lock directories
and atomic os.replace() after verifying exact byte count.
"""
import concurrent.futures
import os
import subprocess
import sys
import threading
import time
import arvados
import arvados.collection
import arvados.keep
from cwl_utils.image_puller import SingularityImagePuller

WORK_ROOT = os.environ.get("NELLY_LOCAL_ROOT", "/ibex/scratch/projects/c2014/nelly-local")
SIF_CACHE = os.path.join(WORK_ROOT, "sif_cache")
tls = threading.local()


def get_kc() -> arvados.keep.KeepClient:
    if not hasattr(tls, "kc"):
        api = arvados.api("v1")
        tls.kc = arvados.keep.KeepClient(api_client=api)
    return tls.kc


def fast_stage(pdh: str, rel_path: str, dest: str, expected_size: int) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest) and os.path.getsize(dest) == expected_size:
        print(f"[fast-stage] Already exists: {dest}", flush=True)
        return
    lock_dir = dest + ".lock"
    try:
        os.mkdir(lock_dir)
    except FileExistsError:
        print(f"[fast-stage] Lock held by another job, skipping: {dest}", flush=True)
        return
    tmp_dest = f"{dest}.fastpart.{os.getpid()}"
    try:
        if os.path.exists(dest) and os.path.getsize(dest) == expected_size:
            return
        t0 = time.time()
        cr = arvados.collection.CollectionReader(pdh)
        fobj = cr[rel_path]
        segs = fobj.segments()
        fd = os.open(tmp_dest, os.O_CREAT | os.O_WRONLY | os.O_TRUNC, 0o644)
        os.ftruncate(fd, expected_size)
        by_loc = {}
        for s in segs:
            by_loc.setdefault(s.locator, []).append((s.range_start, s.segment_offset, s.range_size))
        locs = list(by_loc.keys())
        print(f"[fast-stage] Downloading {dest} ({expected_size} B across {len(locs)} blocks) ...", flush=True)
        done = 0
        lock = threading.Lock()

        def fetch_one(loc: str) -> int:
            nonlocal done
            kc = get_kc()
            blk = memoryview(kc.get(loc))
            for f_off, b_off, length in by_loc[loc]:
                os.pwrite(fd, blk[b_off : b_off + length], f_off)
            with lock:
                done += 1
                if done % 25 == 0 or done == len(locs):
                    elapsed = time.time() - t0
                    mb_s = (done * 64.0) / max(elapsed, 0.1)
                    print(
                        f"[fast-stage] {os.path.basename(dest)}: {done}/{len(locs)} blocks ({elapsed:.1f}s, ~{mb_s:.1f} MiB/s)",
                        flush=True,
                    )
            return len(by_loc[loc])

        with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
            list(ex.map(fetch_one, locs))
        os.close(fd)
        actual = os.path.getsize(tmp_dest)
        if actual != expected_size:
            raise RuntimeError(f"Size mismatch for {dest}: {actual} != {expected_size}")
        os.replace(tmp_dest, dest)
        print(f"[fast-stage] Completed {dest} in {time.time() - t0:.1f}s", flush=True)
    finally:
        if os.path.exists(tmp_dest):
            try:
                os.remove(tmp_dest)
            except OSError:
                pass
        try:
            os.rmdir(lock_dir)
        except OSError:
            pass


def pull_sif(docker_uri: str) -> None:
    os.makedirs(SIF_CACHE, exist_ok=True)
    sif_name = SingularityImagePuller(docker_uri, SIF_CACHE, "singularity", False).get_image_name()
    dst = os.path.join(SIF_CACHE, sif_name)
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        print(f"[fast-stage] SIF already exists: {dst} ({os.path.getsize(dst)} B)", flush=True)
        return
    lock_dir = os.path.join(SIF_CACHE, f".lock-{sif_name}")
    try:
        os.mkdir(lock_dir)
    except FileExistsError:
        print(f"[fast-stage] SIF lock held, skipping: {dst}", flush=True)
        return
    tmp_dst = f"{dst}.tmp.{os.getpid()}"
    try:
        print(f"[fast-stage] Pulling docker://{docker_uri} -> {dst} ...", flush=True)
        subprocess.run(["singularity", "pull", "--force", tmp_dst, f"docker://{docker_uri}"], check=True)
        os.replace(tmp_dst, dst)
        print(f"[fast-stage] Pulled SIF {dst} ({os.path.getsize(dst)} B)", flush=True)
    finally:
        if os.path.exists(tmp_dst):
            try:
                os.remove(tmp_dst)
            except OSError:
                pass
        try:
            os.rmdir(lock_dir)
        except OSError:
            pass


def main() -> int:
    for img in ["quay.io/biocontainers/vg:1.54.0--h9ee0642_0", "google/deepvariant:1.6.0"]:
        pull_sif(img)

    targets = [
        ("1d5ad2aefb6e809dddaffea3a0558415+5657", "ksa_hg38.dist", f"{WORK_ROOT}/inputs/graph/ksa_hg38.dist", 1567026720),
        ("56e94ed8e0d67d0c8e2616959d23ef2e+23423", "ksa_hg38.min", f"{WORK_ROOT}/inputs/graph/ksa_hg38.min", 37348221952),
        ("4bc83e4703035a3eeb142eb17fcc1ce2+44786", "HG003_R2.fastq.gz", f"{WORK_ROOT}/inputs/giab/hg003/HG003_R2.fastq.gz", 33451687186),
        ("0c559a512e87c566c3a3de02d056b4a5+33970", "HG002.novaseq.pcr-free.30x.R2.fastq.gz", f"{WORK_ROOT}/inputs/giab/hg002/HG002.novaseq.pcr-free.30x.R2.fastq.gz", 26634286830),
        ("4bc83e4703035a3eeb142eb17fcc1ce2+44786", "HG003_R1.fastq.gz", f"{WORK_ROOT}/inputs/giab/hg003/HG003_R1.fastq.gz", 31969557485),
        ("0c559a512e87c566c3a3de02d056b4a5+33970", "HG002.novaseq.pcr-free.30x.R1.fastq.gz", f"{WORK_ROOT}/inputs/giab/hg002/HG002.novaseq.pcr-free.30x.R1.fastq.gz", 25664092287),
    ]
    for pdh, rel, dst, sz in targets:
        fast_stage(pdh, rel, dst, sz)
    return 0


if __name__ == "__main__":
    sys.exit(main())
