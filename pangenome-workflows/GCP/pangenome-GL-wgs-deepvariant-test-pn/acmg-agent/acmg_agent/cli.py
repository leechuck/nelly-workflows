"""Command-line entrypoint for acmg-agent (`extract` and `assess` subcommands)."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import sys
import pysam

from acmg_agent.candidates import (
    choose_leading,
    rank_key,
    read_eh_expansions,
    read_tsv,
    select_candidates,
)
from acmg_agent.agent import run_all
from acmg_agent.classify import classify
from acmg_agent.criteria import evaluate_all
from acmg_agent.evidence import build_evidence
from acmg_agent.llm_backends import select_backend_from_env
from acmg_agent.inputs import (
    ClinvarIndex,
    VepIndex,
    load_gene_table,
    norm_chrom,
    read_slivar_info,
)
from acmg_agent.repeats import load_strchive
from acmg_agent.model import Evidence
from acmg_agent.report import (
    build_deterministic_narrative,
    serialize_candidate_assessment,
    write_json,
    write_markdown,
)


def _read_fai_contigs(fai_path: Path | None) -> list[tuple[str, int]]:
    if not fai_path or not fai_path.exists():
        return []
    contigs: list[tuple[str, int]] = []
    for line in fai_path.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) >= 2:
            try:
                contigs.append((cols[0].strip(), int(cols[1].strip())))
            except ValueError:
                continue
    return contigs


def cmd_extract(args: argparse.Namespace) -> int:
    genes = read_tsv(args.genes_tsv) if args.genes_tsv and Path(args.genes_tsv).exists() else []
    variants = read_tsv(args.variants_tsv) if args.variants_tsv and Path(args.variants_tsv).exists() else []
    cands = select_candidates(genes, variants, top=args.top)

    fai_contigs = _read_fai_contigs(Path(args.fai) if args.fai else None)
    fai_names = {c[0] for c in fai_contigs}
    prefer_chr = any(name.startswith("chr") for name in fai_names) if fai_names else True

    header_lines = ["##fileformat=VCFv4.2"]
    for cname, clen in fai_contigs:
        header_lines.append(f"##contig=<ID={cname},length={clen}>")
    header_lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")

    records: list[tuple[str, int, str, str, str, str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    for c in cands:
        raw_chrom = str(c.get("chrom") or "").strip()
        pos_str = str(c.get("pos") or "").strip()
        ref = str(c.get("ref") or "").strip()
        alt = str(c.get("alt") or "").strip()
        if not raw_chrom or not pos_str or not ref or not alt:
            continue
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        bare = norm_chrom(raw_chrom)
        chrom = f"chr{bare}" if (prefer_chr and bare not in ("MT",)) else bare
        if fai_names and chrom not in fai_names and raw_chrom in fai_names:
            chrom = raw_chrom
        key = (chrom, pos, ref, alt)
        if key in seen:
            continue
        seen.add(key)
        rsid = str(c.get("rsid") or ".").strip() or "."
        gt = str(c.get("gt") or "0/1").strip() or "0/1"
        records.append((chrom, pos, rsid, ref, alt, gt))

    records.sort(key=lambda r: (r[0], r[1]))
    out_path = Path(args.output_vcf)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        for hl in header_lines:
            f.write(hl + "\n")
        for chrom, pos, rsid, ref, alt, gt in records:
            f.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    return 0


def _candidates_from_vep_vcf(vep_path: Path, gene_table: dict[str, dict] | None, vep_idx: VepIndex | None) -> list[dict]:
    cands: list[dict] = []
    with pysam.VariantFile(str(vep_path)) as vf:
        for rec in vf:
            chrom = str(rec.chrom)
            pos = int(rec.pos)
            ref = str(rec.ref)
            for alt in (rec.alts or ()):
                alt_s = str(alt)
                annot = vep_idx.lookup(chrom, pos, ref, alt_s) if vep_idx else None
                gene = annot.gene if annot and annot.gene else "UNKNOWN"
                gt_row = gene_table.get(gene) if gene_table and gene else None
                moi = (gt_row.get("moi") if gt_row else None) or "OTHER"
                gt_str = (annot.gt if annot and annot.gt else None) or "0/1"
                hgvs_parts = [p for p in (gene, annot.transcript if annot else None, annot.hgvsc if annot else None, annot.hgvsp if annot else None) if p]
                cands.append({
                    "exomiser_rank": "",
                    "gene": gene,
                    "moi": moi,
                    "combined": "",
                    "pheno_score": "",
                    "human_pheno_score": "",
                    "var_score": "",
                    "chrom": chrom,
                    "pos": str(pos),
                    "ref": ref,
                    "alt": alt_s,
                    "gt": gt_str,
                    "effect": annot.consequence if annot else "",
                    "hgvs": ":".join(hgvs_parts),
                    "rsid": str(rec.id or "") if rec.id != "." else "",
                    "clinvar": "",
                    "max_freq": "",
                    "path_score": "",
                    "disease": "",
                    "pheno_evidence": "",
                })
    return cands


_CLASS_RANK = {
    "Pathogenic": 0,
    "Likely pathogenic": 1,
    "VUS": 2,
    "Likely benign": 3,
    "Benign": 4,
}


def _candidate_key(c: dict) -> tuple[str, str, str, str, str]:
    return (
        str(c.get("gene", "")),
        str(c.get("chrom", "")),
        str(c.get("pos", "")),
        str(c.get("ref", "")),
        str(c.get("alt", "")),
    )


def cmd_assess(args: argparse.Namespace) -> int:
    errors: list[str] = []
    llm_info = {"backend": "none", "model": None, "status": "disabled"}

    try:
        vep_idx = None
        if args.vep_vcf:
            try:
                vep_idx = VepIndex.from_vcf(args.vep_vcf)
            except Exception as exc:
                errors.append(f"Failed to parse VEP VCF ({args.vep_vcf}): {exc}")

        clinvar_idx = None
        if args.clinvar_vcf:
            try:
                clinvar_idx = ClinvarIndex(args.clinvar_vcf)
            except Exception as exc:
                errors.append(f"Failed to open ClinVar VCF ({args.clinvar_vcf}): {exc}")

        gene_table = None
        if args.gene_table:
            try:
                gene_table = load_gene_table(args.gene_table)
            except Exception as exc:
                errors.append(f"Failed to load gene table ({args.gene_table}): {exc}")

        slivar_info = None
        if args.slivar_vcf:
            try:
                slivar_info = read_slivar_info(args.slivar_vcf)
            except Exception as exc:
                errors.append(f"Failed to read slivar VCF ({args.slivar_vcf}): {exc}")

        strchive_loci = None
        strchive_used = False
        if getattr(args, "strchive_json", None):
            try:
                strchive_loci = load_strchive(args.strchive_json)
                strchive_used = True
            except Exception as exc:
                errors.append(f"Failed to load STRchive JSON ({args.strchive_json}): {exc}")

        expansions: list[dict] = []
        if args.eh_json:
            try:
                expansions = read_eh_expansions(args.eh_json, strchive_loci=strchive_loci)
            except Exception as exc:
                errors.append(f"Failed to read ExpansionHunter JSON ({args.eh_json}): {exc}")

        variants_rows: list[dict] = []
        if args.from_vcf or not args.genes_tsv:
            raw_candidates = (
                _candidates_from_vep_vcf(Path(args.vep_vcf), gene_table, vep_idx)
                if args.vep_vcf and Path(args.vep_vcf).exists() and vep_idx is not None
                else []
            )
        else:
            genes_rows = read_tsv(args.genes_tsv)
            variants_rows = read_tsv(args.variants_tsv) if args.variants_tsv and Path(args.variants_tsv).exists() else []
            raw_candidates = select_candidates(genes_rows, variants_rows, top=args.top)

        # Deterministic pass first: the WP4 agent layer adjusts these criteria
        # afterwards under tool-call/time budgets, fail-closed.
        pairs: list[tuple[dict, Evidence, list, int, tuple[str, str, str, str, str]]] = []
        for c in raw_candidates:
            ev = build_evidence(
                c,
                variants_rows=variants_rows,
                vep=vep_idx,
                clinvar=clinvar_idx,
                gene_table=gene_table,
                slivar_info=slivar_info,
            )
            if ev.clinvar_sig and (c.get("clinvar") or "").strip().upper() in ("", "NOT_PROVIDED", "NOT_AVAILABLE", "."):
                c["clinvar"] = ev.clinvar_sig.upper()
            det = evaluate_all(ev)
            if args.from_vcf or not args.genes_tsv:
                tier = _CLASS_RANK.get(classify(det).class_, 2)
            else:
                tier = rank_key(c)[0]
            pairs.append((c, ev, det, tier, _candidate_key(c)))

        # --- WP4 LLM layer: one bounded review pass (tier 0-2 + top 5, <= 8
        # variants, 12 tool calls and 6 min each, 50 min global, contract
        # enforced in agent.validate_submission; rejected/failed/timeout
        # fall back to the deterministic criteria). ------------------------
        llm_status = "disabled"
        assistant: list[dict] | None = None
        review_order: list[int]
        leading_exp_obj: dict | None = None
        if args.from_vcf or not args.genes_tsv:
            review_order = sorted(
                range(len(pairs)),
                key=lambda i: (
                    _CLASS_RANK.get(classify(pairs[i][2]).class_, 9),
                    -classify(pairs[i][2]).points,
                    str(pairs[i][0].get("gene", "")),
                ),
            )
            path_exps = [e for e in expansions if e.get("status") == "pathogenic_range"]
            if path_exps and (not pairs or _CLASS_RANK.get(classify(pairs[review_order[0]][2]).class_, 9) > 0):
                leading_exp_obj = path_exps[0]
        else:
            lead_pre, _ = choose_leading(raw_candidates, expansions)
            if lead_pre is not None and lead_pre[0] == "expansion":
                leading_exp_obj = lead_pre[1]
            lead_snv_key = (
                _candidate_key(lead_pre[1])
                if lead_pre is not None and lead_pre[0] == "snv"
                else None
            )
            review_order = sorted(
                range(len(pairs)),
                key=lambda i: (
                    0 if (lead_snv_key is not None and pairs[i][4] == lead_snv_key) else 1,
                    rank_key(pairs[i][0]),
                    _CLASS_RANK.get(classify(pairs[i][2]).class_, 9),
                    -classify(pairs[i][2]).points,
                ),
            )

        for e in expansions:
            e["agent_status"] = "not_reviewed"
            e["adjusted_by_agent"] = []
            e["agent_rationale"] = None
            e["agent_narrative"] = None

        if not args.no_llm:
            if getattr(args, "sample", None):
                os.environ.setdefault("ACMG_SAMPLE", str(args.sample))
            try:
                backend = select_backend_from_env()
            except Exception as exc:
                backend = None
                llm_status = "failed"
                errors.append(f"LLM backend selection failed: {exc}")
            if backend is not None:
                llm_info["backend"] = str(getattr(backend, "backend_name", "llama"))
                llm_info["model"] = getattr(backend, "model_name", None)
                for e in expansions:
                    if e is not leading_exp_obj:
                        e["agent_status"] = "skipped"
                ordered = [pairs[i] for i in review_order]
                agent_cands = [
                    {
                        "evidence": ev,
                        "criteria": list(det),
                        "tier": tier,
                        "rank": r,
                        "is_leading": (r == 0 and leading_exp_obj is None),
                        "gene": ev.gene,
                    }
                    for r, (_c, ev, det, tier, _k) in enumerate(ordered)
                ]
                try:
                    _reviewed, llm_status = run_all(
                        agent_cands,
                        backend,
                        leading_expansion=leading_exp_obj,
                    )
                    assistant = agent_cands
                except Exception as exc:
                    errors.append(f"LLM layer failed: {exc}")
                    llm_status = "failed"
                llm_info["status"] = llm_status

        # Serialise in review order; reviewed candidates carry the agent's
        # (validated) criteria plus adjusted_by_agent flags.
        assessed_by_key: dict[tuple[str, str, str, str, str], dict] = {}
        assessed_list: list[dict] = []
        for idx_pos, pair_idx in enumerate(review_order):
            c, ev, det, tier, key = pairs[pair_idx]
            if assistant is not None:
                agent = assistant[idx_pos]
                final_criteria = agent["criteria"]
                adjusted = list(agent.get("adjusted") or [])
                agent_status = str(agent.get("agent_status") or "skipped")
                agent_rationale = agent.get("agent_rationale")
                agent_narrative = agent.get("narrative") if agent_status == "ok" else None
            else:
                final_criteria = det
                adjusted = []
                agent_status = "not_reviewed"
                agent_rationale = None
                agent_narrative = None
            clf = classify(final_criteria)
            item = serialize_candidate_assessment(
                c,
                ev,
                final_criteria,
                clf,
                tier=tier,
                adjusted_codes=adjusted,
                agent_rationale=agent_rationale,
                agent_narrative=agent_narrative,
                agent_status=agent_status,
            )
            assessed_list.append(item)
            assessed_by_key[key] = item

        if args.from_vcf or not args.genes_tsv:
            assessed_list.sort(key=lambda a: (_CLASS_RANK.get(a["classification"], 9), -a["points"], a["gene"]))
            if leading_exp_obj is not None:
                exp = leading_exp_obj
                leading_dict = {
                    "kind": "expansion",
                    "gene": exp["gene"],
                    "variant": f"{exp['gene']} repeat ({exp['genotypes']})",
                    "hgvs": None,
                    "zygosity": "het" if exp["allele1"] != exp["allele2"] else "hom",
                    "moi": exp["moi"],
                    "disease": exp["disease"],
                    "classification": "Pathogenic",
                    "points": None,
                    "rule_based_class": "Pathogenic",
                    "conflict": False,
                    "allele1": exp["allele1"],
                    "allele2": exp["allele2"],
                    "genotypes": exp["genotypes"],
                    "ci": exp["ci"],
                    "normal_max": exp["normal_max"],
                    "pathogenic_min": exp["pathogenic_min"],
                    "ranges": exp.get("ranges"),
                    "status": exp["status"],
                    "agent_status": exp.get("agent_status", "not_reviewed"),
                    "adjusted_by_agent": list(exp.get("adjusted_by_agent") or []),
                    "agent_rationale": exp.get("agent_rationale"),
                    "agent_narrative": exp.get("agent_narrative"),
                }
            else:
                leading_dict = assessed_list[0] if assessed_list else None
            near_ties_list: list[dict] = []
        else:
            lead_choice, alt_cands = choose_leading(raw_candidates, expansions)
            # Sort assessed_list by Exomiser/blurb rank_key so leading SNV is first
            order_keys = [
                (c.get("gene", ""), str(c.get("chrom", "")), str(c.get("pos", "")), str(c.get("ref", "")), str(c.get("alt", "")))
                for c in sorted(raw_candidates, key=rank_key)
            ]
            assessed_list = [assessed_by_key[k] for k in order_keys if k in assessed_by_key]

            if lead_choice is None:
                leading_dict = None
            elif lead_choice[0] == "expansion":
                exp = lead_choice[1]
                leading_dict = {
                    "kind": "expansion",
                    "gene": exp["gene"],
                    "variant": f"{exp['gene']} repeat ({exp['genotypes']})",
                    "hgvs": None,
                    "zygosity": "het" if exp["allele1"] != exp["allele2"] else "hom",
                    "moi": exp["moi"],
                    "disease": exp["disease"],
                    "classification": "Pathogenic",
                    "points": None,
                    "rule_based_class": "Pathogenic",
                    "conflict": False,
                    "allele1": exp["allele1"],
                    "allele2": exp["allele2"],
                    "genotypes": exp["genotypes"],
                    "ci": exp["ci"],
                    "normal_max": exp["normal_max"],
                    "pathogenic_min": exp["pathogenic_min"],
                    "ranges": exp.get("ranges"),
                    "status": exp["status"],
                    "agent_status": exp.get("agent_status", "not_reviewed"),
                    "adjusted_by_agent": list(exp.get("adjusted_by_agent") or []),
                    "agent_rationale": exp.get("agent_rationale"),
                    "agent_narrative": exp.get("agent_narrative"),
                }
            else:
                lead_c = lead_choice[1]
                lkey = (
                    lead_c.get("gene", ""),
                    str(lead_c.get("chrom", "")),
                    str(lead_c.get("pos", "")),
                    str(lead_c.get("ref", "")),
                    str(lead_c.get("alt", "")),
                )
                lead_item = assessed_by_key.get(lkey)
                if lead_item is not None:
                    leading_dict = {
                        "kind": "snv",
                        "gene": lead_item["gene"],
                        "variant": lead_item["variant"],
                        "hgvs": lead_item["hgvs"],
                        "zygosity": lead_item["zygosity"],
                        "moi": lead_item["moi"],
                        "disease": lead_item["disease"],
                        "classification": lead_item["classification"],
                        "points": lead_item["points"],
                        "rule_based_class": lead_item["rule_based_class"],
                        "conflict": lead_item["conflict"],
                        "agent_status": lead_item.get("agent_status"),
                        "adjusted_by_agent": lead_item.get("adjusted_by_agent", []),
                        "agent_rationale": lead_item.get("agent_rationale"),
                        "agent_narrative": lead_item.get("agent_narrative"),
                    }
                else:
                    leading_dict = None

            near_ties_list = []
            for ac in alt_cands:
                akey = (
                    ac.get("gene", ""),
                    str(ac.get("chrom", "")),
                    str(ac.get("pos", "")),
                    str(ac.get("ref", "")),
                    str(ac.get("alt", "")),
                )
                if akey in assessed_by_key:
                    aitem = assessed_by_key[akey]
                    near_ties_list.append({
                        "gene": aitem["gene"],
                        "variant": aitem["variant"],
                        "classification": aitem["classification"],
                        "points": aitem["points"],
                    })

        # Populate llm.reviewed, llm.unreviewed, and llm.statuses in rank order (leading first)
        reviewed_genes: list[str] = []
        unreviewed_genes: list[str] = []
        status_map: dict[str, str] = {}
        ordered_items_for_llm: list[dict] = []
        if isinstance(leading_dict, dict) and leading_dict.get("kind") == "expansion":
            ordered_items_for_llm.append(leading_dict)
        ordered_items_for_llm.extend(assessed_list)
        for item in ordered_items_for_llm:
            g = str(item.get("gene") or "")
            st = str(item.get("agent_status") or "not_reviewed")
            if not g:
                continue
            if g not in status_map or st == "ok":
                status_map[g] = st
            if st == "ok":
                if g not in reviewed_genes:
                    reviewed_genes.append(g)
        for item in ordered_items_for_llm:
            g = str(item.get("gene") or "")
            if g and g not in reviewed_genes and g not in unreviewed_genes:
                unreviewed_genes.append(g)
        llm_info["reviewed"] = reviewed_genes
        llm_info["unreviewed"] = unreviewed_genes
        llm_info["statuses"] = status_map

        narrative = build_deterministic_narrative(leading_dict, near_ties_list, assessed_list, expansions)
        top_agent_narrative = (
            leading_dict.get("agent_narrative")
            if isinstance(leading_dict, dict)
            else (assessed_list[0].get("agent_narrative") if assessed_list else None)
        )
        phenotype = (args.phenotype or "")
        data = {
            "schema_version": 1,
            "sample": args.sample,
            "strchive_used": strchive_used,
            "sex": args.sex or None,
            "hpo_terms": [t for t in re.split(r"[\s,]+", args.hpo_terms or "") if t],
            "phenotype": phenotype or None,
            "llm": llm_info,
            "leading_candidate": leading_dict,
            "near_tie_alternatives": near_ties_list,
            "candidates": assessed_list,
            "expansions": expansions,
            "narrative": narrative,
            "agent_narrative": top_agent_narrative,
            "errors": errors,
        }
    except Exception as exc:
        errors.append(f"Unexpected error during assessment: {exc}")
        data = {
            "schema_version": 1,
            "sample": getattr(args, "sample", "proband"),
            "llm": llm_info,
            "leading_candidate": None,
            "near_tie_alternatives": [],
            "candidates": [],
            "expansions": [],
            "narrative": f"Assessment encountered an error: {exc}",
            "errors": errors,
        }

    write_json(args.output_json, data)
    write_markdown(args.output_md, data)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="acmg-agent", description="ACMG/AMP 2015 + ClinGen SVI deterministic engine & CLI")
    sub = p.add_subparsers(dest="subcommand", required=True)

    ext = sub.add_parser("extract", help="Extract top candidate variants from Exomiser TSVs into a mini-VCF for VEP")
    ext.add_argument("--genes-tsv", required=True)
    ext.add_argument("--variants-tsv", required=True)
    ext.add_argument("--fai", default=None)
    ext.add_argument("--top", type=int, default=15)
    ext.add_argument("--output-vcf", required=True)

    ass = sub.add_parser("assess", help="Evaluate ACMG criteria and emit JSON + Markdown reports")
    ass.add_argument("--genes-tsv", default=None)
    ass.add_argument("--variants-tsv", default=None)
    ass.add_argument("--eh-json", default=None)
    ass.add_argument("--strchive-json", default=None)
    ass.add_argument("--vep-vcf", required=True)
    ass.add_argument("--clinvar-vcf", required=True)
    ass.add_argument("--gene-table", required=True)
    ass.add_argument("--slivar-vcf", default=None)
    ass.add_argument("--from-vcf", action="store_true", help="Evaluate all variants in --vep-vcf without Exomiser TSVs")
    ass.add_argument("--no-llm", action="store_true", help="Run deterministic ACMG engine only")
    ass.add_argument("--sample", default="proband")
    ass.add_argument("--hpo-terms", default="", help="HPO term IDs (comma or whitespace separated), free-text phenotype labels go through env ACMG_PHENOTYPE")
    # Phenotype labels arrive via the ACMG_PHENOTYPE env var (CWL
    # EnvVarRequirement), keeping shell/hostile characters out of argv.
    ass.add_argument("--phenotype", default=os.environ.get("ACMG_PHENOTYPE", ""))
    ass.add_argument("--sex", default="")
    ass.add_argument("--top", type=int, default=15)
    ass.add_argument("--output-md", required=True)
    ass.add_argument("--output-json", required=True)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    if args.subcommand == "extract":
        return cmd_extract(args)
    if args.subcommand == "assess":
        return cmd_assess(args)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
