#!/usr/bin/env python3
"""Real-model smoke test for the WP4 LLM layer (run inside nelly-acmg-agent:v1 on ws).

Builds Evidence/Criterion objects from the WP3 TS2 fixtures via
evidence.build_evidence + criteria.evaluate_all (the engine CLI is used the
same way it is wired in production), then runs run_llm_layer (llama-cpp
backend from ACMG_MODEL_PATH / ACMG_THREADS) and prints:
  * the candidate review set and per-candidate llm statuses,
  * the wall time and llm.status,
  * the HEXA criteria table after the agent pass.

Usage (from the repo checkout mounted at /work):
    python /work/acmg-agent/scripts/smoke_llm.py \
        --fixtures /work/acmg-agent/tests/fixtures/acmg \
        --gene-table /work/gene_table.tsv \
        --sample TS2 [--no-llm] [--threads 8] [--top 10]

No network is needed beyond the model file baked into the image.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

_here = Path(__file__).resolve().parent.parent
# in-image at /app, ws smoke run with the repo mounted at /work, or direct
for _p in (Path("/app"), Path("/work/acmg-agent"), _here):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from acmg_agent.agent import run_all  # noqa: E402
from acmg_agent.candidates import rank_key, read_tsv, select_candidates  # noqa: E402
from acmg_agent.criteria import evaluate_all  # noqa: E402
from acmg_agent.evidence import build_evidence  # noqa: E402
from acmg_agent.inputs import ClinvarIndex, VepIndex, load_gene_table  # noqa: E402
from acmg_agent.llm_backends import select_backend_from_env  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="smoke-test run_llm_layer on the TS2 fixtures")
    ap.add_argument("--fixtures", default="/work/acmg-agent/tests/fixtures/acmg")
    ap.add_argument("--gene-table", required=True)
    ap.add_argument("--genes-tsv", default=None, help="default: <fixtures>/TS2/exomiser.genes.tsv")
    ap.add_argument("--variants-tsv", default=None, help="default: <fixtures>/TS2/exomiser.variants.tsv")
    ap.add_argument("--sample", default="TS2")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--lead-gene", default="HEXA", help="gene whose criteria table is printed")
    ap.add_argument("--no-llm", action="store_true")
    ap.add_argument("--max-variants", type=int, default=8)
    ap.add_argument("--global-seconds", type=int, default=3000)
    args = ap.parse_args(argv)

    fix = Path(args.fixtures)
    genes_path = Path(args.genes_tsv or fix / args.sample / "exomiser.genes.tsv")
    variants_path = Path(args.variants_tsv or fix / args.sample / "exomiser.variants.tsv")

    genes_rows = read_tsv(genes_path)
    variants_rows = read_tsv(variants_path)
    raw_candidates = select_candidates(genes_rows, variants_rows, top=args.top)
    vep_vcf = fix / "vep.vcf"

    vep_idx = None
    if vep_vcf.exists():
        vep_idx = VepIndex.from_vcf(vep_vcf)
    clinvar_idx = ClinvarIndex(fix / "clinvar.vcf.gz")
    gene_table = load_gene_table(args.gene_table)

    # rank order as in the workflow (rank_key), Evidence + deterministic criteria
    order = sorted(range(len(raw_candidates)), key=lambda i: rank_key(raw_candidates[i]))
    pack = []
    for rank, i in enumerate(order):
        c = raw_candidates[i]
        ev = build_evidence(c, variants_rows=variants_rows, vep=vep_idx,
                            clinvar=clinvar_idx, gene_table=gene_table)
        det = evaluate_all(ev)
        met = [f"{c2.code}({c2.strength})" for c2 in det if c2.status == "met"]
        print(f"[deterministic] rank {rank + 1}: {ev.gene} {ev.chrom}:{ev.pos} {ev.ref}>{ev.alt} "
              f"{ev.consequence} met={met}")
        pack.append({"evidence": ev, "criteria": list(det), "tier": rank_key(c)[0], "rank": rank})

    if args.no_llm:
        backend_name, llm_status = "none", "off"
        for p in pack:
            p["agent_status"] = "skipped"
    else:
        backend = select_backend_from_env()
        backend_name = getattr(backend, "backend_name", "none")
        model = getattr(backend, "model_name", None)
        llm_status = "off" if backend is None else None
        print(f"[backend] {backend_name} model={model} threads={os.environ.get('ACMG_THREADS', '?')} "
              f"ctx={os.environ.get('ACMG_N_CTX', '8192')}")
        t0 = time.monotonic()
        pack, llm_status = run_all(pack, backend, global_seconds=args.global_seconds,
                                   max_variants=args.max_variants)
        wall = time.monotonic() - t0
        print(f"[wall time] {wall:.1f}s (per-variant budget 360s, global {args.global_seconds}s)")

    statuses = {p["evidence"].gene: p.get("agent_status") for p in pack}
    print(f"[llm.status] {llm_status}")
    print(f"[agent statuses] {statuses}")

    for p in pack:
        ev = p["evidence"]
        if ev.gene != args.lead_gene:
            continue
        print(f"[{args.lead_gene} criteria table] adjusted_by_agent codes: {p.get('adjusted')}")
        print(f"[{args.lead_gene} agent_rationale] {p.get('agent_rationale')}")
        for c in p["criteria"]:
            if c.status in ("met", "not_evaluable"):
                flag = " *" if c.code in (p.get("adjusted") or []) else ""
                print(f"  {c.code:<6} {c.status:<14} {str(c.strength or '-'):<12} {c.rationale}{flag}")

    out = {
        "sample": args.sample,
        "backend": backend_name,
        "llm_status": llm_status,
        "agent_statuses": statuses,
        "candidates": [
            {"gene": p["evidence"].gene, "agent_status": p["agent_status"],
             "adjusted": p.get("adjusted"), "rationale": p.get("agent_rationale")}
            for p in pack
        ],
    }
    print("[json]", json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
