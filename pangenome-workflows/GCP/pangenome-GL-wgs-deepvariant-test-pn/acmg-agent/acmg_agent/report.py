"""JSON (schema_version 1) and Markdown report formatting for ACMG assessments."""
from __future__ import annotations

import json
from pathlib import Path
from acmg_agent.candidates import format_expansions_table
from acmg_agent.model import Classification, Criterion, Evidence


def _format_variant_str(chrom: str, pos: int | str, ref: str, alt: str) -> str:
    c = str(chrom).strip()
    if c and not c.startswith("chr") and c not in ("MT",):
        c = f"chr{c}"
    return f"{c}:{pos} {ref}>{alt}"


def serialize_criterion(c: Criterion) -> dict:
    return {
        "code": c.code,
        "status": c.status,
        "strength": c.strength,
        "rationale": c.rationale,
    }


def serialize_candidate_assessment(
    cand: dict,
    ev: Evidence,
    criteria: list[Criterion],
    clf: Classification,
    tier: int,
    llm_adjustments: list[dict] | None = None,
    rationale: str | None = None,
) -> dict:
    met_codes = [
        f"{c.code} ({c.strength})" if c.strength else c.code
        for c in criteria
        if c.status == "met"
    ]
    default_rationale = (
        f"{ev.gene} ({ev.moi}) {_format_variant_str(ev.chrom, ev.pos, ev.ref, ev.alt)} "
        f"({ev.zygosity}) classified as {clf.class_} ({clf.points} pts; rule-based: {clf.rule_based_class}). "
        f"Met criteria: {', '.join(met_codes) if met_codes else 'none'}."
    )
    return {
        "kind": "snv",
        "gene": ev.gene,
        "variant": _format_variant_str(ev.chrom, ev.pos, ev.ref, ev.alt),
        "hgvs": cand.get("hgvs") or ev.hgvsc or ev.hgvsp or "",
        "zygosity": ev.zygosity,
        "moi": ev.moi,
        "disease": ev.disease,
        "tier": tier,
        "classification": clf.class_,
        "points": clf.points,
        "rule_based_class": clf.rule_based_class,
        "conflict": clf.conflict,
        "evidence": {
            "transcript": ev.transcript,
            "consequence": ev.consequence,
            "mane_select": ev.mane_select,
            "lof_mechanism": ev.lof_mechanism,
            "max_af": ev.max_af,
            "popmax_af": ev.popmax_af,
            "nhomalt": ev.nhomalt,
            "revel": ev.revel,
            "spliceai": ev.spliceai,
            "clinvar_id": ev.clinvar_id,
            "clinvar_sig": ev.clinvar_sig,
            "clinvar_stars": ev.clinvar_stars,
            "exomiser_pheno_score": ev.human_pheno_score,
        },
        "criteria": [serialize_criterion(c) for c in criteria],
        "llm_adjustments": llm_adjustments or [],
        "rationale": rationale or default_rationale,
    }


def build_deterministic_narrative(
    leading: dict | None,
    near_ties: list[dict],
    candidates: list[dict],
    expansions: list[dict],
) -> str:
    if leading is None:
        return (
            "No confident pathogenic or likely pathogenic candidate was identified among the "
            "evaluated SNV/indel variants or short tandem repeat loci."
        )
    if leading.get("kind") == "expansion":
        return (
            f"Leading finding: pathogenic-range repeat expansion in **{leading['gene']}** "
            f"(alleles `{leading.get('genotypes', '')}`, 95% CI `{leading.get('ci', '')}`; "
            f"pathogenic threshold >= {leading.get('pathogenic_min')} repeats; MOI {leading.get('moi')}) "
            f"associated with **{leading.get('disease')}**."
        )
    lead_str = (
        f"Leading candidate: **{leading['gene']}** `{leading['variant']}` (`{leading.get('hgvs', '')}`, "
        f"zygosity `{leading['zygosity']}`, MOI `{leading['moi']}`), associated with **{leading.get('disease') or 'unspecified disorder'}**. "
        f"ACMG/AMP classification: **{leading['classification']}** ({leading['points']} points; "
        f"rule-based: {leading['rule_based_class']})."
    )
    if near_ties:
        alt_genes = ", ".join(f"**{a['gene']}** ({a.get('classification', 'VUS')})" for a in near_ties)
        lead_str += f" Near-tie differential candidates requiring clinical correlation: {alt_genes}."
    return lead_str


def render_markdown(data: dict) -> str:
    lines: list[str] = ["# ACMG/AMP Variant Interpretation Report", ""]
    if data.get("strchive_used") is False:
        lines.append(
            "_Note: `--strchive-json` was not provided; repeat expansion thresholds fell back to the built-in legacy table._"
        )
        lines.append("")

    narrative = data.get("narrative") or ""
    lines.append("## 1. Leading Conclusion")
    lines.append("")
    lines.append(narrative)
    lines.append("")

    lead = data.get("leading_candidate")
    cands = data.get("candidates") or []

    # Full 28-criteria table for the leading SNV candidate (or top SNV candidate if expansion leads)
    featured_snv = None
    if lead and lead.get("kind") == "snv":
        featured_snv = next((c for c in cands if c["gene"] == lead["gene"] and c["variant"] == lead["variant"]), None)
    if featured_snv is None and cands:
        featured_snv = cands[0]

    if featured_snv is not None:
        lines.append(
            f"## 2. ACMG 28-Criteria Evaluation — {featured_snv['gene']} ({featured_snv['variant']})"
        )
        lines.append("")
        lines.append(
            f"- **Classification**: **{featured_snv['classification']}** "
            f"({featured_snv['points']} pts; rule-based: {featured_snv['rule_based_class']}; "
            f"conflict: `{featured_snv['conflict']}`)"
        )
        lines.append(f"- **HGVS**: `{featured_snv.get('hgvs') or '-'}`")
        lines.append(f"- **Zygosity / MOI**: `{featured_snv['zygosity']}` / `{featured_snv['moi']}`")
        lines.append(f"- **Disease**: {featured_snv.get('disease') or '-'}")
        lines.append("")
        lines.append("| Criterion | Status | Strength | Rationale |")
        lines.append("|---|---|---|---|")
        for c in featured_snv.get("criteria", []):
            lines.append(
                f"| {c['code']} | {c['status']} | {c['strength'] or '—'} | {c['rationale']} |"
            )
        lines.append("")

    lines.append("## 3. Assessed SNV / Indel Candidates")
    lines.append("")
    if cands:
        lines.append("| # | Gene | Variant | Zygosity | MOI | Class | Points | Met Criteria | Disease |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for idx, c in enumerate(cands, 1):
            met = [
                f"{cr['code']}({cr['strength']})" if cr.get("strength") else cr["code"]
                for cr in c.get("criteria", [])
                if cr.get("status") == "met"
            ]
            lines.append(
                f"| {idx} | {c['gene']} | {c['variant']} | {c['zygosity']} | {c['moi']} | "
                f"**{c['classification']}** | {c['points']} | {', '.join(met) or '—'} | "
                f"{(c.get('disease') or '-')[:45]} |"
            )
    else:
        lines.append("_No SNV/indel candidates evaluated._")
    lines.append("")

    lines.append("## 4. Short Tandem Repeat Expansions (ExpansionHunter)")
    lines.append("")
    lines.append(format_expansions_table(data.get("expansions") or []))
    lines.append("")

    if data.get("errors"):
        lines.append("## 5. Warnings / Input Notes")
        lines.append("")
        for err in data["errors"]:
            lines.append(f"- {err}")
        lines.append("")

    return "\n".join(lines)


def write_json(path: Path | str, data: dict) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2) + "\n")


def write_markdown(path: Path | str, data: dict) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(render_markdown(data) + "\n")
