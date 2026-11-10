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


def serialize_criterion(c: Criterion, adjusted_by_agent: bool = False) -> dict:
    return {
        "code": c.code,
        "status": c.status,
        "strength": c.strength,
        "rationale": c.rationale,
        "adjusted_by_agent": adjusted_by_agent,
    }


def serialize_candidate_assessment(
    cand: dict,
    ev: Evidence,
    criteria: list[Criterion],
    clf: Classification,
    tier: int,
    llm_adjustments: list[dict] | None = None,
    rationale: str | None = None,
    adjusted_codes: list[str] | None = None,
    agent_rationale: str | None = None,
    agent_narrative: str | None = None,
    agent_status: str | None = None,
) -> dict:
    if adjusted_codes is not None:
        adjusted_list = [str(code) for code in adjusted_codes]
    elif llm_adjustments is not None:
        adjusted_list = [
            str(a["code"])
            for a in llm_adjustments
            if isinstance(a, dict) and a.get("code")
        ]
    else:
        adjusted_list = []
    adjusted_set = set(adjusted_list)
    computed_adjustments = (
        llm_adjustments
        if llm_adjustments is not None
        else [
            {
                "code": c.code,
                "status": c.status,
                "strength": c.strength,
                "reason": c.rationale,
            }
            for c in criteria
            if c.code in adjusted_set
        ]
    )
    if agent_status:
        resolved_status = str(agent_status)
    elif agent_rationale or agent_narrative or adjusted_list or llm_adjustments:
        resolved_status = "ok"
    else:
        resolved_status = "not_reviewed"
    met_codes = [
        f"{c.code} ({c.strength})" if c.strength else c.code
        for c in criteria
        if c.status == "met"
    ]
    founder_note = ""
    if ev.max_af_any_pop is not None and (ev.max_af is None or ev.max_af_any_pop > ev.max_af):
        pop = ev.max_af_any_pop_pop or ev.max_af_any_pop_source or "founder"
        founder_note = f" (founder-population AF {ev.max_af_any_pop:.4f} ({pop}))"
    default_rationale = (
        f"{ev.gene} ({ev.moi}) {_format_variant_str(ev.chrom, ev.pos, ev.ref, ev.alt)} "
        f"({ev.zygosity}) classified as {clf.class_} ({clf.points} pts; rule-based: {clf.rule_based_class}). "
        f"Met criteria: {', '.join(met_codes) if met_codes else 'none'}{founder_note}."
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
            "max_af_source": ev.max_af_source,
            "max_af_any_pop": ev.max_af_any_pop,
            "max_af_any_pop_source": ev.max_af_any_pop_source,
            "max_af_any_pop_pop": ev.max_af_any_pop_pop,
            "popmax_af": ev.popmax_af,
            "nhomalt": ev.nhomalt,
            "revel": ev.revel,
            "spliceai": ev.spliceai,
            "clinvar_id": ev.clinvar_id,
            "clinvar_sig": ev.clinvar_sig,
            "clinvar_stars": ev.clinvar_stars,
            "exomiser_pheno_score": ev.human_pheno_score,
        },
        "criteria": [
            serialize_criterion(c, adjusted_by_agent=(c.code in adjusted_set))
            for c in criteria
        ],
        "llm_adjustments": computed_adjustments,
        "adjusted_by_agent": adjusted_list,
        "agent_status": resolved_status,
        "agent_rationale": agent_rationale,
        "agent_narrative": agent_narrative,
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

    llm = data.get("llm") or {}
    llm_backend = llm.get("backend") or "none"
    llm_model = llm.get("model") or "none"
    llm_status = llm.get("status") or "disabled"

    lead = data.get("leading_candidate")
    cands = data.get("candidates") or []

    featured_snv = None
    if lead and lead.get("kind") == "snv":
        featured_snv = next((c for c in cands if c["gene"] == lead["gene"] and c["variant"] == lead["variant"]), None)
    if featured_snv is None and cands:
        featured_snv = cands[0]

    lead_agent_status = (
        (lead.get("agent_status") if isinstance(lead, dict) else None)
        or (featured_snv.get("agent_status") if featured_snv else None)
    )
    lead_gene = (
        (lead.get("gene") if isinstance(lead, dict) else None)
        or (featured_snv.get("gene") if featured_snv else None)
        or "leading candidate"
    )

    lines.append(
        f"- **LLM Review**: backend=`{llm_backend}`, model=`{llm_model}`, status=`{llm_status}`"
    )
    if llm_status in ("disabled", "off"):
        lines.append(
            f"> **NOTE — Deterministic-only report (`llm.status={llm_status}`)**: "
            f"Leading candidate (**{lead_gene}**) was NOT reviewed by an LLM."
        )
    elif lead_agent_status != "ok" or llm_status in ("partial", "timeout", "failed", "rejected"):
        if lead_agent_status != "ok":
            lines.append(
                f"> **WARNING — Leading candidate NOT reviewed by LLM "
                f"(llm.status=`{llm_status}`, leading candidate agent_status=`{lead_agent_status or 'unreviewed'}`)**: "
                f"The leading candidate (**{lead_gene}**) did not complete LLM review; its classification and "
                f"criteria below are **deterministic-only**."
            )
        else:
            lines.append(
                f"> **NOTE — Partial LLM review (`llm.status={llm_status}`)**: "
                f"Leading candidate (**{lead_gene}**) was reviewed (`agent_status=ok`), "
                f"but one or more secondary candidates were not reviewed."
            )
    lines.append("")

    narrative = data.get("narrative") or ""
    lines.append("## 1. Leading Conclusion")
    lines.append("")
    lines.append(narrative)
    lines.append("")

    lead_narrative = (
        data.get("agent_narrative")
        or (lead.get("agent_narrative") if isinstance(lead, dict) else None)
        or (featured_snv.get("agent_narrative") if (lead is None and featured_snv) else None)
    )
    lead_rationale = (
        (lead.get("agent_rationale") if isinstance(lead, dict) else None)
        or (featured_snv.get("agent_rationale") if (lead is None and featured_snv) else None)
    )
    if lead_narrative:
        lines.append(f"**Leading Candidate LLM Narrative ({lead_gene})**: {lead_narrative}")
        lines.append("")
    elif lead_agent_status == "ok" and lead_rationale:
        lines.append(f"**Leading Candidate LLM Rationale ({lead_gene})**: {lead_rationale}")
        lines.append("")
    elif llm_status not in ("disabled", "off") and lead_agent_status != "ok":
        lines.append(
            f"_Leading candidate **{lead_gene}** was NOT reviewed by the LLM "
            f"(status: `{lead_agent_status or llm_status}`); deterministic summary shown above._"
        )
        lines.append("")

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
        lines.append(f"- **LLM Review Status**: `{featured_snv.get('agent_status') or 'not_reviewed'}`")
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

    lines.append(f"## 5. LLM review ({llm_model}, {llm_status})")
    lines.append("")
    lines.append(
        f"- **LLM configuration (model, status)**: backend=`{llm_backend}`, "
        f"model=`{llm_model}`, status=`{llm_status}`"
    )
    if "reviewed" in llm or "unreviewed" in llm:
        lines.append(f"- **Reviewed genes**: {', '.join(llm.get('reviewed') or []) or 'none'}")
        lines.append(f"- **Unreviewed genes**: {', '.join(llm.get('unreviewed') or []) or 'none'}")
    lines.append("")
    has_lead_exp = isinstance(lead, dict) and lead.get("kind") == "expansion"
    if not cands and not has_lead_exp:
        lines.append("_No candidates evaluated._")
        lines.append("")
    else:
        if has_lead_exp:
            exp_status = lead.get("agent_status") or "not_reviewed"
            lines.append(
                f"### {lead['gene']} (repeat expansion {lead.get('genotypes', '')}) — status: `{exp_status}`"
            )
            lines.append(f"- **Status**: `{exp_status}`")
            lines.append("- **Adjusted criteria**: none (repeat expansion; no ACMG criteria adjustments)")
            lines.append(f"- **Rationale**: {lead.get('agent_rationale') or '—'}")
            lines.append(f"- **Narrative**: {lead.get('agent_narrative') or '—'}")
            lines.append("")
        for c in cands:
            c_status = c.get("agent_status") or "not_reviewed"
            adj_items = c.get("llm_adjustments") or [
                {
                    "code": cr["code"],
                    "status": cr["status"],
                    "strength": cr.get("strength"),
                    "reason": cr.get("rationale", ""),
                }
                for cr in c.get("criteria", [])
                if cr.get("adjusted_by_agent")
            ]
            lines.append(f"### {c['gene']} ({c['variant']}) — status: `{c_status}`")
            lines.append(f"- **Status**: `{c_status}`")
            if adj_items:
                adj_strs = [
                    f"`{a['code']}` -> `{a['status']}` (`{a.get('strength') or '—'}`): {a.get('reason') or a.get('rationale') or ''}"
                    for a in adj_items
                ]
                lines.append(f"- **Adjusted criteria**: {'; '.join(adj_strs)}")
            else:
                lines.append("- **Adjusted criteria**: none")
            lines.append(f"- **Rationale**: {c.get('agent_rationale') or '—'}")
            lines.append(f"- **Narrative**: {c.get('agent_narrative') or '—'}")
            lines.append("")

    if data.get("errors"):
        lines.append("## 6. Warnings / Input Notes")
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
