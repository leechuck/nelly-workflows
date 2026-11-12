"""Clinical PDF report renderer (HTML + CSS -> PDF via WeasyPrint, design §11 + R9 + R10).

Main body presents the deterministic ACMG/AMP 2015 + ClinGen SVI evaluation:
  - Header block (unnumbered: patient/sample ID, sex, HPO terms + phenotype labels,
    GRCh38, pipeline/image versions from env, report date)
  - 1. Clinical Conclusion (leading finding, classification + points + combining
    rule, supported vs inconclusive decision, structured summary, near-tie differential)
  - 2. Primary Variant Evaluation:
       * SNV/indel-led: ACMG/AMP 28-Criteria Evaluation Table for the leading finding
       * Repeat-expansion-led: Primary Repeat Expansion Decision Table (locus, gene,
         alleles + 95% CI, STRchive normal/intermediate/pathogenic ranges + catalogue
         version, decision rule applied, MOI, disease, and explicit statement that
         ACMG/AMP sequence-variant criteria do not apply to repeat expansions),
         followed by 2.1 Highest-ranked SNV/indel candidate (secondary finding —
         not the reported diagnosis)
  - 3. Ranked SNV / Indel Candidates Table (fixed layout with colgroup + border-box,
    fitting A4 portrait with full right margin) + Compact P/LP Table for secondary
    P/LP candidates labelled 'secondary finding — not the reported diagnosis'
  - 4. Short Tandem Repeat Expansions (pathogenic/intermediate first, normal loci
    compact)
  - Appendix A. Automated Reviewer Commentary (summary status table first, followed
    by per-candidate rationale and narrative — explicitly labelled as explanatory
    commentary, not the basis of the clinical decision)
  - Appendix B. Methods, Data Versions, Limitations & Run Identifiers
  - Appendix C. Full 28-Criteria Tables for Secondary Pathogenic / Likely
    Pathogenic Candidates ('secondary finding — not the reported diagnosis')
"""
from __future__ import annotations

from datetime import datetime, timezone
import html
import os
from pathlib import Path
import re
import sys
from typing import Any
from acmg_agent.policy import REDACTION_NOTE

PATHOGENIC_CODES = (
    "PVS1",
    "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
)
BENIGN_CODES = (
    "BA1",
    "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
)

STRENGTH_POINTS = {
    ("P", "very_strong"): 8,
    ("P", "strong"): 4,
    ("P", "moderate"): 2,
    ("P", "supporting"): 1,
    ("B", "stand_alone"): -8,
    ("B", "strong"): -4,
    ("B", "moderate"): -2,
    ("B", "supporting"): -1,
}

KV_COLGROUP = (
    "<colgroup>"
    '<col style="width: 18%"/>'
    '<col style="width: 32%"/>'
    '<col style="width: 18%"/>'
    '<col style="width: 32%"/>'
    "</colgroup>"
)

CRITERIA_COLGROUP = (
    "<colgroup>"
    '<col style="width: 7%"/>'
    '<col style="width: 21%"/>'
    '<col style="width: 15%"/>'
    '<col style="width: 37%"/>'
    '<col style="width: 20%"/>'
    "</colgroup>"
)


def _esc(val: Any, default: str = "—") -> str:
    if val is None or val == "":
        return default
    return html.escape(str(val))


def _md_inline_to_html(val: Any, default: str = "—") -> str:
    """Escape HTML and convert inline Markdown (**bold**, `code`) to HTML tags."""
    if val is None or val == "":
        return default
    escaped = html.escape(str(val))
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    escaped = re.sub(r"`([^`]+)`", r"<code>\1</code>", escaped)
    return escaped


def _format_non_bottlenecked_af(ev: dict) -> str:
    """Format non-bottlenecked AF: numeric only when > 0, else 'not observed in gnomAD v4 (exomes/genomes)'."""
    af_val = ev.get("max_af")
    if isinstance(af_val, (int, float)) and af_val > 0.0:
        af_str = f"{af_val:.2e} ({ev.get('max_af_source') or 'gnomAD v4'})"
    else:
        af_str = "not observed in gnomAD v4 (exomes/genomes)"
    any_pop = ev.get("max_af_any_pop")
    if isinstance(any_pop, (int, float)) and any_pop > 0.0 and (af_val is None or any_pop > af_val):
        af_str += f" [founder pop {ev.get('max_af_any_pop_pop') or 'max'}: {any_pop:.4f}]"
    return af_str


def _criterion_source(code: str, ev: dict) -> str:
    code_u = (code or "").upper()
    if code_u == "PVS1":
        return "VEP 116 CSQ (MANE) + gene_table (LOEUF/LoF)"
    if code_u in ("PS1", "PM5"):
        return "ClinVar GRCh38 VCF + VEP 116 CSQ"
    if code_u in ("PM2", "BA1", "BS1"):
        af_val = ev.get("max_af")
        src = ev.get("max_af_source")
        if isinstance(af_val, (int, float)) and af_val > 0.0 and src:
            return f"gnomAD v4 ({src})"
        return "gnomAD v4 (exomes/genomes)"
    if code_u == "BS2":
        return "slivar gnomAD v4 (gnomad_nhomalt)"
    if code_u == "PM1":
        return "gene_table (hotspot / ClinVar missense)"
    if code_u == "PM3":
        return "VCF genotype + gene_table MOI + gnomAD AF"
    if code_u in ("PM4", "BP3", "BP7"):
        return "VEP 116 CSQ + RepeatMasker / SpliceAI"
    if code_u in ("PP2", "BP1"):
        return "VEP 116 CSQ + gene_table (ClinVar missense/LoF)"
    if code_u in ("PP3", "BP4"):
        return "VEP 116 / Exomiser (REVEL, SpliceAI)"
    if code_u == "PP4":
        return "Exomiser 2512 HPO phenotype match"
    if code_u in ("PP5", "BP6"):
        return "ClinGen SVI (Biesecker et al. 2018)"
    if code_u in ("PS2", "PM6", "PP1", "BS4"):
        return "Single-sample workflow (no pedigree/trio)"
    if code_u in ("PS3", "BS3"):
        return "No curated functional assay DB wired"
    if code_u == "PS4":
        return "ClinVar VCF + gnomAD AF"
    if code_u in ("BP2", "BP5"):
        return "Exomiser candidate set + gene_table MOI"
    return "ACMG deterministic engine"


def _criterion_points(code: str, strength: str | None) -> int:
    is_benign = (code or "").upper().startswith("B")
    direction = "B" if is_benign else "P"
    return STRENGTH_POINTS.get((direction, strength or ""), 0)


def _decision_label(status: str, strength: str | None, code: str) -> tuple[str, str, str]:
    """Return (decision_badge_html, applied_strength_str, css_row_class)."""
    st = (status or "").lower()
    if st == "met":
        pts = _criterion_points(code, strength)
        pts_str = f"+{pts}" if pts > 0 else str(pts)
        str_label = (
            f"{strength or 'met'} ({pts_str} pt{'s' if abs(pts) != 1 else ''})"
            if strength != "stand_alone"
            else "stand_alone"
        )
        badge = '<span class="badge badge-met">Decision made: Met</span>'
        return badge, _esc(str_label), "row-met"
    if st == "not_met":
        badge = '<span class="badge badge-not-met">Decision made: Not met</span>'
        return badge, "— (0 pts)", "row-not-met"
    if st == "not_applicable":
        badge = '<span class="badge badge-na">Not made: Not applicable</span>'
        return badge, "— (retired)", "row-na"
    badge = '<span class="badge badge-ne">Not made: Not evaluable</span>'
    return badge, "—", "row-ne"


def _render_criteria_group_rows(criteria_by_code: dict[str, dict], codes: tuple[str, ...], ev: dict) -> str:
    rows: list[str] = []
    for code in codes:
        c = criteria_by_code.get(code) or {
            "code": code,
            "status": "not_evaluable",
            "strength": None,
            "rationale": "not evaluated",
        }
        badge_html, strength_str, row_cls = _decision_label(
            c.get("status", "not_evaluable"),
            c.get("strength"),
            code,
        )
        rat = _esc(c.get("rationale") or "—")
        src = _esc(_criterion_source(code, ev))
        rows.append(
            f'<tr class="{row_cls}">'
            f'<td class="col-code"><strong>{_esc(code)}</strong></td>'
            f'<td class="col-decision">{badge_html}</td>'
            f'<td class="col-strength">{strength_str}</td>'
            f'<td class="col-evidence">{rat}</td>'
            f'<td class="col-source">{src}</td>'
            f"</tr>"
        )
    return "\n".join(rows)


def _render_candidate_criteria_section(
    cand: dict,
    is_leading: bool = False,
    role_override: str | None = None,
) -> str:
    gene = _esc(cand.get("gene"))
    variant = _esc(cand.get("variant"))
    hgvs = _esc(cand.get("hgvs"))
    zyg = _esc(cand.get("zygosity"))
    moi = _esc(cand.get("moi"))
    disease = _esc(cand.get("disease"))
    clf = _esc(cand.get("classification"))
    pts = cand.get("points")
    pts_str = f"{pts:+d}" if isinstance(pts, int) else _esc(pts, "0")
    rb_clf = _esc(cand.get("rule_based_class"))
    conflict = bool(cand.get("conflict"))
    ev = cand.get("evidence") or {}

    criteria_list = cand.get("criteria") or []
    by_code = {str(c.get("code")): c for c in criteria_list if isinstance(c, dict)}

    path_rows = _render_criteria_group_rows(by_code, PATHOGENIC_CODES, ev)
    benign_rows = _render_criteria_group_rows(by_code, BENIGN_CODES, ev)

    if role_override:
        role_tag = _esc(role_override)
    elif is_leading:
        role_tag = "Leading SNV/Indel Candidate"
    else:
        role_tag = "Secondary Pathogenic / Likely Pathogenic Candidate (secondary finding — not the reported diagnosis)"

    secondary_meta = (
        ""
        if is_leading
        else '<span><strong>Finding Role:</strong> <span class="sec-note">secondary finding — not the reported diagnosis</span></span>'
    )
    conflict_note = (
        ' <span class="badge badge-warn">Pathogenic/Benign Evidence Conflict</span>'
        if conflict
        else ""
    )

    af_str = _format_non_bottlenecked_af(ev)

    return f"""
    <div class="criteria-card">
      <div class="cand-header-group">
        <h3>{role_tag}: {gene} — <code>{variant}</code></h3>
        <div class="cand-meta">
          {secondary_meta}
          <span><strong>HGVS:</strong> <code>{hgvs}</code></span>
          <span><strong>Zygosity / MOI:</strong> <code>{zyg}</code> / <code>{moi}</code></span>
          <span><strong>Disease:</strong> {disease}</span>
          <span><strong>Non-bottlenecked AF:</strong> <code>{_esc(af_str)}</code></span>
          <span><strong>Point Tally:</strong> <strong>{_esc(pts_str)} pts</strong></span>
          <span><strong>Richards Combining Rule:</strong> <strong>{rb_clf}</strong></span>
          <span><strong>Final Classification:</strong> <strong class="clf-label">{clf}</strong>{conflict_note}</span>
        </div>
        <h4>Pathogenic Criteria (PVS1, PS1–PS4, PM1–PM6, PP1–PP5)</h4>
      </div>
      <table class="criteria-table">
        {CRITERIA_COLGROUP}
        <thead>
          <tr>
            <th>Code</th>
            <th>Decision Status</th>
            <th>Applied Strength</th>
            <th>Evidence Value &amp; Evaluation Rationale</th>
            <th>Evidence Source</th>
          </tr>
        </thead>
        <tbody>
          {path_rows}
        </tbody>
      </table>

      <h4>Benign Criteria (BA1, BS1–BS4, BP1–BP7)</h4>
      <table class="criteria-table">
        {CRITERIA_COLGROUP}
        <thead>
          <tr>
            <th>Code</th>
            <th>Decision Status</th>
            <th>Applied Strength</th>
            <th>Evidence Value &amp; Evaluation Rationale</th>
            <th>Evidence Source</th>
          </tr>
        </thead>
        <tbody>
          {benign_rows}
        </tbody>
      </table>
    </div>
    """


def _render_expansion_primary_section(lead: dict, data: dict, lead_snv: dict | None) -> str:
    """Render Section 2 when the leading finding is a short tandem repeat expansion."""
    gene = _esc(lead.get("gene"))
    locus_id = _esc(lead.get("locus_id") or lead.get("gene"))
    genotypes = str(lead.get("genotypes") or "")
    ci = _esc(lead.get("ci"))
    a1 = lead.get("allele1")
    a2 = lead.get("allele2")
    alleles_num = [x for x in (a1, a2) if isinstance(x, int)]
    if not alleles_num and genotypes:
        alleles_num = [int(x) for x in genotypes.split("/") if x.strip().isdigit()]
    max_al = lead.get("max_allele")
    if not isinstance(max_al, int):
        max_al = max(alleles_num) if alleles_num else "—"
    norm_max = lead.get("normal_max")
    path_min = lead.get("pathogenic_min")
    ranges = lead.get("ranges") or {}
    inter_tuple = ranges.get("intermediate")
    if isinstance(inter_tuple, (list, tuple)) and len(inter_tuple) == 2 and inter_tuple[0] is not None and inter_tuple[1] is not None:
        inter_str = f"{inter_tuple[0]}–{inter_tuple[1]} repeats"
    elif isinstance(norm_max, int) and isinstance(path_min, int) and path_min > norm_max + 1:
        inter_str = f"{norm_max + 1}–{path_min - 1} repeats"
    else:
        inter_str = "none (adjacent normal/pathogenic thresholds)"

    strchive_used = bool(data.get("strchive_used"))
    cat_ver = (
        "STRchive v2025.1 locus definitions + ExpansionHunter v5 repeat catalogue (strchive_used=True)"
        if strchive_used
        else "STRchive / ExpansionHunter v5 curated repeat locus catalogue (strchive_used=False)"
    )
    rule_applied = f"larger allele {max_al} \u2265 pathogenic_min {path_min} \u2192 pathogenic range"

    sec2_html = f"""
    <div class="cand-header-group">
      <h3>Leading Finding (Primary Repeat Expansion Decision): {gene} (Locus <code>{locus_id}</code>)</h3>
      <table class="kv-table">
        {KV_COLGROUP}
        <tr>
          <th>Repeat Locus / Gene</th>
          <td>Locus <code>{locus_id}</code> / Gene <strong>{gene}</strong></td>
          <th>Repeat Alleles &amp; 95% CI</th>
          <td>Alleles <code>{_esc(genotypes)}</code> repeats (95% CI: <code>{ci}</code>)</td>
        </tr>
        <tr>
          <th>STRchive Normal Range</th>
          <td>&le; <code>{_esc(norm_max)}</code> repeats</td>
          <th>STRchive Intermediate Range</th>
          <td><code>{_esc(inter_str)}</code></td>
        </tr>
        <tr>
          <th>STRchive Pathogenic Range</th>
          <td>&ge; <code>{_esc(path_min)}</code> repeats</td>
          <th>Reference Catalogue Version</th>
          <td>{_esc(cat_ver)}</td>
        </tr>
        <tr>
          <th>Zygosity / MOI</th>
          <td><code>{_esc(lead.get("zygosity"))}</code> / <code>{_esc(lead.get("moi"))}</code></td>
          <th>Associated Disease</th>
          <td>{_esc(lead.get("disease"))}</td>
        </tr>
        <tr>
          <th>Decision Rule Applied</th>
          <td colspan="3"><strong><code>{_esc(rule_applied)}</code></strong></td>
        </tr>
      </table>
      <div class="notice-box">
        <strong>ACMG/AMP Sequence-Variant Criteria Applicability Notice:</strong>
        The 28 ACMG/AMP 2015 sequence-variant criteria (Richards et al. 2015) are defined for single-nucleotide variants and small indels and <strong>do not apply to short tandem repeat expansions</strong>. Repeat expansions are interpreted against locus-specific repeat-count thresholds, intermediate/reduced-penetrance ranges, and sizing confidence intervals per the ACMG points-to-consider and technical standards for repeat expansion testing (Monaghan et al. 2020, <em>Genet Med</em> 22:981; Bean et al. 2020, <em>Genet Med</em> 22:443).
      </div>
    </div>
    """
    if lead_snv is not None:
        sec2_html += f"""
        <h3>2.1 Highest-ranked SNV/indel candidate (secondary finding)</h3>
        <p><em>Secondary finding — not the reported diagnosis: The 28-criteria ACMG/AMP sequence-variant evaluation below is for the highest-ranked SNV/indel candidate in this sample and is a secondary finding, not the primary reported diagnosis.</em></p>
        {_render_candidate_criteria_section(
            lead_snv,
            is_leading=False,
            role_override="Highest-ranked SNV/indel candidate (secondary finding — not the reported diagnosis)",
        )}
        """
    return sec2_html


def _narrative_is_complete(text: object, status: object) -> bool:
    if str(status or "") == "ok_truncated":
        return False
    if not isinstance(text, str) or not text.strip():
        return False
    return bool(re.search(r'[.!?]["\'”’)\]]*\s*$', text.strip()))



def _redaction_html(item: dict | None) -> str:
    """One-line note when policy removed recommendation sentences from LLM text."""
    if isinstance(item, dict) and item.get("narrative_redacted"):
        return f' <em>({_esc(REDACTION_NOTE)})</em>'
    return ""

def _build_structured_summary_html(data: dict, lead: dict | None, near_ties: list[dict], cands: list[dict]) -> str:
    """Render the Section 1 Summary from structured fields (or convert inline Markdown if complete agent narrative)."""
    if isinstance(lead, dict) and _narrative_is_complete(lead.get("agent_narrative"), lead.get("agent_status")):
        return _md_inline_to_html(lead["agent_narrative"]) + _redaction_html(lead)
    truncated_notice = ""
    if isinstance(lead, dict) and (
        lead.get("agent_status") == "ok_truncated"
        or (lead.get("agent_narrative") and not _narrative_is_complete(lead.get("agent_narrative"), lead.get("agent_status")))
    ):
        truncated_notice = (
            " <em>[Note: Automated reviewer narrative was truncated at the token limit after retry "
            "(<code>agent_status=ok_truncated</code>); deterministic clinical summary shown above and "
            "complete reviewer rationale in Appendix A.]</em>"
        )
    if isinstance(lead, dict) and lead.get("kind") == "expansion":
        return (
            f"Leading candidate is a pathogenic-range short tandem repeat expansion in "
            f"<strong>{_esc(lead.get('gene'))}</strong> (alleles <code>{_esc(lead.get('genotypes'))}</code>, "
            f"95% CI <code>{_esc(lead.get('ci'))}</code>; normal &le; <code>{_esc(lead.get('normal_max'))}</code>, "
            f"pathogenic &ge; <code>{_esc(lead.get('pathogenic_min'))}</code>) associated with "
            f"{_esc(lead.get('disease'))} (<code>{_esc(lead.get('moi'))}</code>)."
            f"{truncated_notice}"
        )
    if isinstance(lead, dict):
        lead_snv = next(
            (c for c in cands if c.get("gene") == lead.get("gene") and c.get("variant") == lead.get("variant")),
            None,
        )
        met_parts: list[str] = []
        if lead_snv:
            for cr in (lead_snv.get("criteria") or []):
                if isinstance(cr, dict) and cr.get("status") == "met":
                    code = _esc(cr.get("code"))
                    st = _esc(cr.get("strength") or "")
                    met_parts.append(f"<code>{code}</code> ({st})" if st != "—" else f"<code>{code}</code>")
        met_str = ", ".join(met_parts) if met_parts else "none"
        pts = lead.get("points")
        pts_str = f"{pts:+d} points" if isinstance(pts, int) else _esc(pts, "0 points")
        summary = (
            f"Leading candidate is <strong>{_esc(lead.get('gene'))}</strong> "
            f"<code>{_esc(lead.get('variant'))}</code> (<code>{_esc(lead.get('hgvs'))}</code>, "
            f"<code>{_esc(lead.get('zygosity'))}</code>, MOI <code>{_esc(lead.get('moi'))}</code>), "
            f"classified as <strong>{_esc(lead.get('classification'))}</strong> "
            f"({_esc(pts_str)}; combining rule: {_esc(lead.get('rule_based_class'))}; "
            f"met criteria: {met_str}) for {_esc(lead.get('disease'))}."
        )
        if near_ties:
            alt_str = ", ".join(
                f"<strong>{_esc(t.get('gene'))}</strong> <code>{_esc(t.get('variant'))}</code> "
                f"({_esc(t.get('classification'))}, {_esc(t.get('points'))} pts)"
                for t in near_ties
            )
            summary += f" Near-tie alternatives within tie margin: {alt_str}."
        return summary + truncated_notice
    return _md_inline_to_html(data.get("narrative") or "No pathogenic or likely pathogenic candidate was identified.")


def _render_compact_plp_table(secondary_plp_cands: list[dict]) -> str:
    """Render compact summary table for non-leading P/LP candidates in Section 3."""
    rows: list[str] = []
    for c in secondary_plp_cands:
        met_items: list[str] = []
        ne_items: list[str] = []
        for cr in (c.get("criteria") or []):
            if not isinstance(cr, dict):
                continue
            code = str(cr.get("code") or "")
            status = str(cr.get("status") or "")
            strength = cr.get("strength")
            if status == "met":
                pts = _criterion_points(code, strength)
                pts_str = f"{pts:+d}" if pts != 0 else "0"
                met_items.append(f"<code>{_esc(code)}</code> ({_esc(strength)}, {pts_str} pt{'s' if abs(pts) != 1 else ''})")
            elif status == "not_evaluable":
                ne_items.append(f"<code>{_esc(code)}</code>")
        pts = c.get("points")
        pts_str = f"{pts:+d} pts" if isinstance(pts, int) else _esc(pts, "0 pts")
        rows.append(
            f'<tr class="row-met">'
            f"<td><strong>{_esc(c.get('gene'))}</strong><br/><code>{_esc(c.get('variant'))}</code><br/><code>{_esc(c.get('hgvs'))}</code></td>"
            f"<td><code>{_esc(c.get('zygosity'))}</code> / <code>{_esc(c.get('moi'))}</code><br/>{_esc(c.get('disease'))}</td>"
            f"<td>{'<br/>'.join(met_items) if met_items else '—'}</td>"
            f"<td>{', '.join(ne_items) if ne_items else '—'}</td>"
            f"<td><strong>{_esc(c.get('classification'))}</strong> ({_esc(pts_str)}; rule: {_esc(c.get('rule_based_class'))})"
            f'<br/><span class="sec-note">secondary finding — not the reported diagnosis</span></td>'
            f"</tr>"
        )
    return f"""
    <h3>Secondary Pathogenic / Likely Pathogenic SNV/Indel Candidates (secondary finding — not the reported diagnosis; Full 28-Row Tables in Appendix C)</h3>
    <table class="data-table">
      <colgroup>
        <col style="width: 23%"/>
        <col style="width: 18%"/>
        <col style="width: 24%"/>
        <col style="width: 15%"/>
        <col style="width: 20%"/>
      </colgroup>
      <thead>
        <tr>
          <th>Gene &amp; Variant (HGVS)</th>
          <th>Zyg / MOI / Disease</th>
          <th>Met Criteria (Strength &amp; Points)</th>
          <th>Not-Evaluable Criteria</th>
          <th>Final Class &amp; Role</th>
        </tr>
      </thead>
      <tbody>
        {"".join(rows)}
      </tbody>
    </table>
    """


def render_pdf_html(data: dict) -> str:
    """Render the clinical report HTML+CSS document from assessment JSON data."""
    raw_sample = str(data.get("sample") or "proband").strip() or "proband"
    raw_patient_id = str(data.get("patient_id") or "").strip() or raw_sample
    sample = _esc(raw_sample)
    patient_id = _esc(raw_patient_id)
    sex = _esc(data.get("sex") or "unspecified")
    hpo_list = data.get("hpo_terms") or []
    hpo_str = ", ".join(str(t) for t in hpo_list) if hpo_list else "none provided"
    phenotype_labels = _esc(data.get("phenotype") or "—")
    report_date = _esc(data.get("report_date") or datetime.now(timezone.utc).strftime("%Y-%m-%d"))

    pipeline_version = _esc(
        data.get("pipeline_version")
        or os.environ.get("ACMG_PIPELINE_VERSION")
        or "nelly-workflows modernise-2026-11"
    )
    image_tag = _esc(
        data.get("image_tag")
        or os.environ.get("ACMG_IMAGE_TAG")
        or "nelly-acmg-agent:v7"
    )

    lead = data.get("leading_candidate")
    near_ties = data.get("near_tie_alternatives") or []
    cands = data.get("candidates") or []
    expansions = data.get("expansions") or []
    llm = data.get("llm") or {}
    llm_backend = _esc(llm.get("backend") or "none")
    llm_model = _esc(llm.get("model") or "none")
    llm_status = _esc(llm.get("status") or "disabled")
    llm_error = _esc(llm.get("error") or "")
    llm_error_html = (
        f'<p class="llm-error"><strong>LLM review not run:</strong> {llm_error}. '
        f'All results in this report are the deterministic engine result.</p>'
        if llm_error else ""
    )

    is_expansion_led = isinstance(lead, dict) and lead.get("kind") == "expansion"

    # Determine whether conclusion is supported or inconclusive
    is_supported = False
    if isinstance(lead, dict):
        if is_expansion_led and lead.get("status") == "pathogenic_range":
            is_supported = True
        elif lead.get("classification") in ("Pathogenic", "Likely pathogenic"):
            is_supported = True
    conclusion_status = "SUPPORTED" if is_supported else "INCONCLUSIVE"
    conclusion_badge_cls = "badge-met" if is_supported else "badge-warn"

    # Section 1: Conclusion block
    if is_expansion_led and isinstance(lead, dict):
        lead_summary_html = f"""
        <table class="kv-table">
          {KV_COLGROUP}
          <tr>
            <th>Leading Finding</th>
            <td><strong>{_esc(lead.get("gene"))}</strong> short tandem repeat expansion (alleles <code>{_esc(lead.get("genotypes"))}</code>, 95% CI <code>{_esc(lead.get("ci"))}</code>)</td>
            <th>Conclusion Status</th>
            <td><span class="badge {conclusion_badge_cls}">{conclusion_status}</span></td>
          </tr>
          <tr>
            <th>HGVS / Repeat Thresholds</th>
            <td>Normal &le; <code>{_esc(lead.get("normal_max"))}</code> repeats; Pathogenic &ge; <code>{_esc(lead.get("pathogenic_min"))}</code> repeats (STRchive / ExpansionHunter)</td>
            <th>Classification</th>
            <td><strong>{_esc(lead.get("classification", "Pathogenic"))}</strong> (Locus threshold rule: {_esc(lead.get("rule_based_class", "Pathogenic"))})</td>
          </tr>
          <tr>
            <th>Zygosity / MOI</th>
            <td><code>{_esc(lead.get("zygosity"))}</code> / <code>{_esc(lead.get("moi"))}</code></td>
            <th>Associated Disease</th>
            <td>{_esc(lead.get("disease"))}</td>
          </tr>
        </table>
        """
    elif isinstance(lead, dict):
        pts = lead.get("points")
        pts_str = f"{pts:+d} points" if isinstance(pts, int) else _esc(pts, "0 points")
        lead_summary_html = f"""
        <table class="kv-table">
          {KV_COLGROUP}
          <tr>
            <th>Leading Finding</th>
            <td><strong>{_esc(lead.get("gene"))}</strong> — <code>{_esc(lead.get("variant"))}</code></td>
            <th>Conclusion Status</th>
            <td><span class="badge {conclusion_badge_cls}">{conclusion_status}</span></td>
          </tr>
          <tr>
            <th>HGVS (c. / p.)</th>
            <td><code>{_esc(lead.get("hgvs"))}</code></td>
            <th>Classification</th>
            <td><strong>{_esc(lead.get("classification"))}</strong> ({_esc(pts_str)}; combining rule: {_esc(lead.get("rule_based_class"))})</td>
          </tr>
          <tr>
            <th>Zygosity / MOI</th>
            <td><code>{_esc(lead.get("zygosity"))}</code> / <code>{_esc(lead.get("moi"))}</code></td>
            <th>Associated Disease</th>
            <td>{_esc(lead.get("disease"))}</td>
          </tr>
        </table>
        """
    else:
        lead_summary_html = (
            '<p><strong>Conclusion Status:</strong> <span class="badge badge-warn">INCONCLUSIVE</span> — '
            "No pathogenic or likely pathogenic candidate was identified.</p>"
        )

    summary_html = _build_structured_summary_html(data, lead, near_ties, cands)

    if near_ties:
        ties_str = ", ".join(
            f"<strong>{_esc(t.get('gene'))}</strong> (<code>{_esc(t.get('variant'))}</code>, {_esc(t.get('classification'))}, {_esc(t.get('points'))} pts)"
            for t in near_ties
        )
        diff_html = f"<p><strong>Near-tie differential candidates:</strong> {ties_str}</p>"
    else:
        diff_html = "<p><strong>Near-tie differential candidates:</strong> none within tie margin.</p>"

    # Identify leading SNV (if SNV-led) or highest-ranked SNV (if expansion-led)
    lead_snv = None
    if isinstance(lead, dict) and lead.get("kind") == "snv":
        lead_snv = next(
            (c for c in cands if c.get("gene") == lead.get("gene") and c.get("variant") == lead.get("variant")),
            None,
        )
    if lead_snv is None and cands:
        lead_snv = cands[0]

    leading_snv_key: tuple[str, str] | None = (
        (str(lead_snv.get("gene")), str(lead_snv.get("variant")))
        if (lead_snv is not None and not is_expansion_led)
        else None
    )

    # Section 2: Primary evaluation
    shown_in_sec2_keys: set[tuple[str, str]] = set()
    if is_expansion_led and isinstance(lead, dict):
        sec2_heading = "2. Primary Variant Evaluation — Short Tandem Repeat Expansion (Leading Finding)"
        sec2_intro = ""
        criteria_html = _render_expansion_primary_section(lead, data, lead_snv)
        if lead_snv is not None:
            shown_in_sec2_keys.add((str(lead_snv.get("gene")), str(lead_snv.get("variant"))))
    else:
        sec2_heading = "2. ACMG/AMP 28-Criteria Evaluation (Leading Finding)"
        sec2_intro = (
            "<p>Each of the 28 ACMG/AMP 2015 criteria is evaluated deterministically with ClinGen SVI calibrations "
            "for the leading finding. Every row states whether a decision was made (<strong>Met</strong> with applied "
            "strength, or <strong>Not met</strong>) or not made (<strong>Not evaluable</strong> with reason, or "
            "<strong>Not applicable</strong> for retired criteria PP5/BP6), together with the underlying evidence value, "
            "source, Tavtigian 2020 point tally, and Richards 2015 combining-rule classification.</p>"
        )
        if lead_snv is not None:
            criteria_html = _render_candidate_criteria_section(lead_snv, is_leading=True)
            shown_in_sec2_keys.add((str(lead_snv.get("gene")), str(lead_snv.get("variant"))))
        else:
            criteria_html = "<p><em>No SNV/indel candidates were evaluated against the 28 ACMG/AMP criteria.</em></p>"

    # All non-leading P/LP SNV/indel candidates (labelled 'secondary finding — not the reported diagnosis')
    secondary_plp_cands: list[dict] = []
    appendix_c_cands: list[dict] = []
    for c in cands:
        key = (str(c.get("gene")), str(c.get("variant")))
        if leading_snv_key is not None and key == leading_snv_key:
            continue
        if c.get("classification") in ("Pathogenic", "Likely pathogenic"):
            secondary_plp_cands.append(c)
            if key not in shown_in_sec2_keys:
                appendix_c_cands.append(c)

    # Section 3: Ranked SNV / Indel Candidates table (fixed layout + colgroup fitting A4 portrait) + compact secondary P/LP table
    if cands:
        cand_rows: list[str] = []
        for idx, c in enumerate(cands, 1):
            key = (str(c.get("gene")), str(c.get("variant")))
            is_c_leading = (leading_snv_key is not None and key == leading_snv_key)
            is_secondary_plp = (
                not is_c_leading
                and c.get("classification") in ("Pathogenic", "Likely pathogenic")
            )
            sec_tag_html = (
                '<br/><span class="sec-note">secondary finding — not the reported diagnosis</span>'
                if is_secondary_plp
                else ""
            )
            met_codes = [
                f"<code>{_esc(cr['code'])}</code> ({_esc(cr['strength'])})"
                if cr.get("strength")
                else f"<code>{_esc(cr['code'])}</code>"
                for cr in (c.get("criteria") or [])
                if isinstance(cr, dict) and cr.get("status") == "met"
            ]
            met_cell_html = "<br/>".join(met_codes) if met_codes else "—"
            cand_rows.append(
                f"<tr>"
                f"<td>{idx}</td>"
                f"<td><strong>{_esc(c.get('gene'))}</strong></td>"
                f"<td><code>{_esc(c.get('variant'))}</code></td>"
                f"<td><code>{_esc(c.get('hgvs'))}</code></td>"
                f"<td><code>{_esc(c.get('zygosity'))}</code></td>"
                f"<td><code>{_esc(c.get('moi'))}</code></td>"
                f"<td><strong>{_esc(c.get('classification'))}</strong>{sec_tag_html}</td>"
                f"<td>{_esc(c.get('points'))}</td>"
                f"<td>{met_cell_html}</td>"
                f"<td>{_esc(c.get('disease'))}</td>"
                f"</tr>"
            )
        ranked_table_html = f"""
        <table class="data-table ranked-table">
          <colgroup>
            <col style="width: 3.5%"/>
            <col style="width: 8.5%"/>
            <col style="width: 13%"/>
            <col style="width: 21%"/>
            <col style="width: 4.5%"/>
            <col style="width: 4.5%"/>
            <col style="width: 13.5%"/>
            <col style="width: 4%"/>
            <col style="width: 12.5%"/>
            <col style="width: 15%"/>
          </colgroup>
          <thead>
            <tr>
              <th>#</th>
              <th>Gene</th>
              <th>Variant</th>
              <th>HGVS</th>
              <th>Zyg</th>
              <th>MOI</th>
              <th>Classification</th>
              <th>Pts</th>
              <th>Met Criteria</th>
              <th>Associated Disease</th>
            </tr>
          </thead>
          <tbody>
            {"".join(cand_rows)}
          </tbody>
        </table>
        """
        if secondary_plp_cands:
            ranked_table_html += "\n" + _render_compact_plp_table(secondary_plp_cands)
    else:
        ranked_table_html = "<p><em>No SNV/indel candidates evaluated.</em></p>"

    # Section 4: Repeat expansions (pathogenic/intermediate first, normal compact)
    abnormal_exp = [e for e in expansions if e.get("status") in ("pathogenic_range", "intermediate_or_carrier")]
    normal_exp = [e for e in expansions if e.get("status") not in ("pathogenic_range", "intermediate_or_carrier")]

    exp_parts: list[str] = []
    if abnormal_exp:
        erows: list[str] = []
        for e in abnormal_exp:
            erows.append(
                f'<tr class="row-met">'
                f"<td><strong>{_esc(e.get('gene'))}</strong></td>"
                f"<td><code>{_esc(e.get('genotypes'))}</code></td>"
                f"<td><code>{_esc(e.get('ci'))}</code></td>"
                f"<td>&le; {_esc(e.get('normal_max'))}</td>"
                f"<td>&ge; {_esc(e.get('pathogenic_min'))}</td>"
                f"<td><code>{_esc(e.get('moi'))}</code></td>"
                f"<td><strong>{_esc(e.get('status'))}</strong></td>"
                f"<td>{_esc(e.get('disease'))}</td>"
                f"</tr>"
            )
        exp_parts.append(f"""
        <h4>Pathogenic &amp; Intermediate Repeat Expansion Calls</h4>
        <table class="data-table">
          <colgroup>
            <col style="width: 9%"/>
            <col style="width: 12%"/>
            <col style="width: 14%"/>
            <col style="width: 9%"/>
            <col style="width: 10%"/>
            <col style="width: 7%"/>
            <col style="width: 15%"/>
            <col style="width: 24%"/>
          </colgroup>
          <thead>
            <tr>
              <th>Gene</th>
              <th>Repeat Alleles</th>
              <th>95% CI</th>
              <th>Normal</th>
              <th>Pathogenic</th>
              <th>MOI</th>
              <th>Status</th>
              <th>Associated Disease</th>
            </tr>
          </thead>
          <tbody>
            {"".join(erows)}
          </tbody>
        </table>
        """)
    else:
        exp_parts.append("<p><em>No pathogenic-range or intermediate repeat expansions detected.</em></p>")

    if normal_exp:
        compact_items = [
            f"<strong>{_esc(e.get('gene'))}</strong>: <code>{_esc(e.get('genotypes'))}</code> (normal &le; {_esc(e.get('normal_max'))})"
            for e in normal_exp
        ]
        exp_parts.append(
            f"<p class=\"compact-loci\"><strong>Normal-range STR loci ({len(normal_exp)} loci):</strong> "
            + "; ".join(compact_items)
            + "</p>"
        )
    expansions_html = "\n".join(exp_parts)

    # Appendix A: Automated reviewer commentary — FIRST a summary table (gene, status, adjusted criteria, wall time), THEN rationale/narrative
    reviewed_genes = ", ".join(llm.get("reviewed") or []) or "none"
    unreviewed_genes = ", ".join(llm.get("unreviewed") or []) or "none"
    appendix_a_summary_rows: list[str] = []
    appendix_a_cards: list[str] = []

    if is_expansion_led and isinstance(lead, dict):
        exp_st = _esc(lead.get("agent_status") or "not_reviewed")
        exp_wt = lead.get("agent_wall_seconds")
        exp_wt_str = f"{exp_wt:.2f}s" if isinstance(exp_wt, (int, float)) else "—"
        appendix_a_summary_rows.append(
            f"<tr>"
            f"<td><strong>{_esc(lead.get('gene'))}</strong></td>"
            f"<td>Repeat Expansion <code>{_esc(lead.get('genotypes'))}</code></td>"
            f"<td><code>{exp_st}</code></td>"
            f"<td>none (locus threshold)</td>"
            f"<td><code>{_esc(exp_wt_str)}</code></td>"
            f"</tr>"
        )
        if lead.get("agent_rationale") or lead.get("agent_narrative") or lead.get("agent_status") not in (None, "not_reviewed", "skipped"):
            appendix_a_cards.append(f"""
            <div class="llm-card">
              <h4>{_esc(lead.get("gene"))} (Repeat Expansion <code>{_esc(lead.get("genotypes"))}</code>) — Status: <code>{exp_st}</code></h4>
              <ul>
                <li><strong>Reviewer rationale:</strong> {_md_inline_to_html(lead.get("agent_rationale"))}</li>
                <li><strong>Reviewer narrative:</strong> {_md_inline_to_html(lead.get("agent_narrative"))}{_redaction_html(lead)}</li>
              </ul>
            </div>
            """)

    for c in cands:
        cst = _esc(c.get("agent_status") or "not_reviewed")
        wt = c.get("agent_wall_seconds")
        wt_str = f"{wt:.2f}s" if isinstance(wt, (int, float)) else "—"
        adj_items = c.get("llm_adjustments") or [
            {
                "code": cr.get("code"),
                "status": cr.get("status"),
                "strength": cr.get("strength"),
                "reason": cr.get("rationale"),
            }
            for cr in (c.get("criteria") or [])
            if isinstance(cr, dict) and cr.get("adjusted_by_agent")
        ]
        if adj_items:
            adj_str = "; ".join(
                f"<code>{_esc(a.get('code'))}</code> &rarr; <code>{_esc(a.get('status'))}</code> "
                f"(<code>{_esc(a.get('strength'))}</code>): {_esc(a.get('reason') or a.get('rationale'))}"
                for a in adj_items
            )
            adj_short = ", ".join(
                f"<code>{_esc(a.get('code'))}</code> ({_esc(a.get('strength'))})"
                for a in adj_items
            )
        else:
            adj_str = "none"
            adj_short = "none"
        appendix_a_summary_rows.append(
            f"<tr>"
            f"<td><strong>{_esc(c.get('gene'))}</strong></td>"
            f"<td><code>{_esc(c.get('variant'))}</code></td>"
            f"<td><code>{cst}</code></td>"
            f"<td>{adj_short}</td>"
            f"<td><code>{_esc(wt_str)}</code></td>"
            f"</tr>"
        )
        has_card_content = bool(
            adj_items or c.get("agent_rationale") or c.get("agent_narrative")
            or c.get("agent_status") not in (None, "not_reviewed", "skipped")
        )
        if has_card_content:
            raw_nar = c.get("agent_narrative")
            if c.get("agent_status") == "ok_truncated" or (raw_nar and not _narrative_is_complete(raw_nar, c.get("agent_status"))):
                nar_html = (
                    "<em>[Truncated at token limit after retry (<code>status=ok_truncated</code>); "
                    "refer to complete reviewer rationale above.]</em> "
                    + _md_inline_to_html(raw_nar)
                )
            else:
                nar_html = _md_inline_to_html(raw_nar)
            appendix_a_cards.append(f"""
            <div class="llm-card">
              <h4>{_esc(c.get("gene"))} (<code>{_esc(c.get("variant"))}</code>) — Status: <code>{cst}</code></h4>
              <ul>
                <li><strong>Adjusted criteria:</strong> {adj_str}</li>
                <li><strong>Reviewer rationale:</strong> {_md_inline_to_html(c.get("agent_rationale"))}</li>
                <li><strong>Reviewer narrative:</strong> {nar_html}{_redaction_html(c)}</li>
              </ul>
            </div>
            """)

    if appendix_a_summary_rows:
        appendix_a_table = f"""
        <table class="data-table">
          <colgroup>
            <col style="width: 14%"/>
            <col style="width: 28%"/>
            <col style="width: 16%"/>
            <col style="width: 28%"/>
            <col style="width: 14%"/>
          </colgroup>
          <thead>
            <tr>
              <th>Gene</th>
              <th>Variant / Repeat Locus</th>
              <th>Review Status</th>
              <th>Adjusted Criteria</th>
              <th>Wall Time</th>
            </tr>
          </thead>
          <tbody>
            {"".join(appendix_a_summary_rows)}
          </tbody>
        </table>
        """
        if appendix_a_cards:
            appendix_a_body = appendix_a_table + "\n" + "\n".join(appendix_a_cards)
        else:
            appendix_a_body = (
                appendix_a_table
                + "\n<p><em>No per-candidate LLM rationale or narrative text was generated (see Review Status column above).</em></p>"
            )
    else:
        appendix_a_body = "<p><em>No candidates evaluated.</em></p>"

    # Appendix C: Full 28-criteria tables for additional secondary P/LP candidates
    if appendix_c_cands:
        app_c_sections = "\n".join(
            _render_candidate_criteria_section(c, is_leading=False) for c in appendix_c_cands
        )
        appendix_c_html = f"""
        <div class="appendix-break"></div>
        <h2>Appendix C. Full 28-Criteria Tables for Secondary Pathogenic / Likely Pathogenic Candidates (secondary finding — not the reported diagnosis)</h2>
        <p>Full 28-row ACMG/AMP 2015 + ClinGen SVI evaluations for additional Pathogenic and Likely pathogenic candidates summarized in Section 3. Each variant below is a <strong>secondary finding — not the reported diagnosis</strong>.</p>
        {app_c_sections}
        """
    else:
        appendix_c_html = ""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>ACMG/AMP Clinical Variant Report — Patient ID: {patient_id} (VCF sample: {sample})</title>
  <style>
    @page {{
      size: A4;
      margin: 15mm 12mm 16mm 12mm;
      @bottom-left {{
        content: "ACMG/AMP Clinical Variant Interpretation Report — Patient ID: {patient_id} (VCF sample: {sample}) (GRCh38)";
        font-family: "DejaVu Sans", sans-serif;
        font-size: 7.5pt;
        color: #555555;
      }}
      @bottom-right {{
        content: "Page " counter(page) " of " counter(pages);
        font-family: "DejaVu Sans", sans-serif;
        font-size: 7.5pt;
        color: #555555;
      }}
    }}
    * {{
      box-sizing: border-box;
    }}
    body {{
      font-family: "DejaVu Sans", sans-serif;
      font-size: 8.5pt;
      line-height: 1.35;
      color: #1a1a1a;
      margin: 0;
    }}
    code {{
      font-family: "DejaVu Sans Mono", monospace;
      font-size: 7.8pt;
      overflow-wrap: break-word;
      word-wrap: break-word;
    }}
    h1 {{
      font-size: 14pt;
      margin: 0 0 6px 0;
      color: #0f2c59;
      border-bottom: 2px solid #0f2c59;
      padding-bottom: 4px;
    }}
    h2 {{
      font-size: 11pt;
      margin: 12px 0 6px 0;
      color: #0f2c59;
      border-bottom: 1px solid #cbd5e1;
      padding-bottom: 2px;
      page-break-after: avoid;
      break-after: avoid;
    }}
    h3 {{
      font-size: 9.5pt;
      margin: 8px 0 4px 0;
      color: #1e293b;
      page-break-after: avoid;
      break-after: avoid;
    }}
    h4 {{
      font-size: 8.5pt;
      margin: 6px 0 3px 0;
      color: #334155;
      page-break-after: avoid;
      break-after: avoid;
    }}
    p {{
      margin: 4px 0;
    }}
    .kv-table, .data-table, .criteria-table {{
      width: 100%;
      table-layout: fixed;
      border-collapse: collapse;
      margin: 4px 0 8px 0;
      font-size: 7.8pt;
    }}
    .ranked-table {{
      font-size: 7.4pt;
      line-height: 1.25;
    }}
    .ranked-table code {{
      font-size: 7pt;
    }}
    tr {{
      page-break-inside: avoid;
      break-inside: avoid;
    }}
    .kv-table th, .kv-table td,
    .data-table th, .data-table td,
    .criteria-table th, .criteria-table td {{
      border: 1px solid #cbd5e1;
      padding: 3px 4px;
      text-align: left;
      vertical-align: top;
      overflow-wrap: break-word;
      word-wrap: break-word;
    }}
    .kv-table th {{
      background-color: #f1f5f9;
      font-weight: bold;
    }}
    .data-table th, .criteria-table th {{
      background-color: #1e3a5f;
      color: #ffffff;
      font-weight: bold;
    }}
    .row-met {{
      background-color: #ecfdf5;
    }}
    .row-not-met {{
      background-color: #ffffff;
    }}
    .row-ne {{
      background-color: #f8fafc;
      color: #475569;
    }}
    .row-na {{
      background-color: #f1f5f9;
      color: #64748b;
    }}
    .badge {{
      display: inline-block;
      padding: 1px 5px;
      border-radius: 3px;
      font-size: 7.2pt;
      font-weight: bold;
    }}
    .badge-met {{
      background-color: #065f46;
      color: #ffffff;
    }}
    .badge-not-met {{
      background-color: #e2e8f0;
      color: #1e293b;
    }}
    .badge-ne {{
      background-color: #fef3c7;
      color: #92400e;
    }}
    .badge-na {{
      background-color: #e2e8f0;
      color: #64748b;
    }}
    .badge-warn {{
      background-color: #991b1b;
      color: #ffffff;
    }}
    .sec-note {{
      color: #92400e;
      font-size: 7pt;
      font-style: italic;
      font-weight: normal;
    }}
    .cand-header-group {{
      page-break-inside: avoid;
      break-inside: avoid;
      page-break-after: avoid;
      break-after: avoid;
    }}
    .cand-meta {{
      background-color: #f8fafc;
      border: 1px solid #cbd5e1;
      padding: 4px 6px;
      margin-bottom: 6px;
      font-size: 8pt;
      page-break-inside: avoid;
      break-inside: avoid;
      page-break-after: avoid;
      break-after: avoid;
    }}
    .cand-meta span {{
      display: inline-block;
      margin-right: 12px;
    }}
    .criteria-card {{
      margin-bottom: 10px;
    }}
    .compact-loci {{
      font-size: 7.5pt;
      color: #334155;
      background-color: #f8fafc;
      border: 1px solid #e2e8f0;
      padding: 4px 6px;
    }}
    .appendix-break {{
      page-break-before: always;
    }}
    .llm-error {{ color: #8a1c1c; border-left: 3px solid #8a1c1c; padding-left: 6px; }}
    .notice-box {{
      background-color: #fffbeb;
      border: 1px solid #f59e0b;
      padding: 6px 8px;
      margin: 6px 0;
      font-size: 8pt;
    }}
    .llm-card {{
      border: 1px solid #cbd5e1;
      padding: 4px 8px;
      margin-bottom: 6px;
      background-color: #f8fafc;
      page-break-inside: avoid;
      break-inside: avoid;
    }}
    .llm-card ul {{
      margin: 3px 0;
      padding-left: 16px;
    }}
  </style>
</head>
<body>

  <h1>Clinical Variant Interpretation Report (ACMG/AMP 2015 + ClinGen SVI)</h1>
  <table class="kv-table">
    {KV_COLGROUP}
    <tr>
      <th>Patient / Sample ID</th>
      <td>Patient ID: <code>{patient_id}</code> (VCF sample: <code>{sample}</code>)</td>
      <th>Sex (Metadata)</th>
      <td><code>{sex}</code></td>
    </tr>
    <tr>
      <th>HPO Terms</th>
      <td><code>{_esc(hpo_str)}</code></td>
      <th>Phenotype Labels</th>
      <td>{phenotype_labels}</td>
    </tr>
    <tr>
      <th>Reference Build</th>
      <td><code>GRCh38 (hg38)</code></td>
      <th>Report Date (UTC)</th>
      <td><code>{report_date}</code></td>
    </tr>
    <tr>
      <th>Pipeline / Image</th>
      <td colspan="3"><code>{pipeline_version}</code> (step image: <code>{image_tag}</code>)</td>
    </tr>
  </table>

  <h2>1. Clinical Conclusion</h2>
  {lead_summary_html}
  <p><strong>Summary:</strong> {summary_html}</p>
  {diff_html}

  <h2>{sec2_heading}</h2>
  {sec2_intro}
  {criteria_html}

  <h2>3. Ranked SNV / Indel Candidates</h2>
  {ranked_table_html}

  <h2>4. Short Tandem Repeat Expansions (ExpansionHunter + STRchive)</h2>
  {expansions_html}

  <div class="appendix-break"></div>
  <h2>Appendix A. Automated Reviewer Commentary ({llm_model}, {llm_status})</h2>
  <div class="notice-box">
    <strong>Explanatory Commentary Notice:</strong> This appendix contains the automated reviewer commentary (backend: <code>{llm_backend}</code>, model: <code>{llm_model}</code>, run status: <code>{llm_status}</code>). It is provided strictly as an explanation and audit trail, <strong>not as the basis of the clinical classification decision</strong>. All classifications in the main report are computed deterministically from the 28 ACMG/AMP criteria under fail-closed evidence-pack validation rules.
  </div>
  {llm_error_html}
  <p><strong>Reviewed genes:</strong> <code>{_esc(reviewed_genes)}</code> &nbsp;|&nbsp; <strong>Unreviewed genes:</strong> <code>{_esc(unreviewed_genes)}</code></p>
  {appendix_a_body}

  <h2>Appendix B. Methods, Data Versions, Limitations &amp; Run Identifiers</h2>
  <table class="kv-table">
    {KV_COLGROUP}
    <tr>
      <th>ACMG/AMP Framework</th>
      <td colspan="3">Richards et al. 2015 (<em>Genet Med</em> 17:405) 28-criteria combining rules + Tavtigian et al. 2020 Bayesian point scale (PVS=8, PS=4, PM=2, PP=1, BS=&minus;4, BP=&minus;1; Pathogenic &ge;10 pts, Likely pathogenic 6–9 pts, VUS 0–5 pts, Likely benign &minus;1 to &minus;5 pts, Benign &le;&minus;6 pts or BA1). Short tandem repeat expansions are evaluated by locus-specific thresholds per ACMG repeat-expansion points-to-consider (Monaghan et al. 2020; Bean et al. 2020).</td>
    </tr>
    <tr>
      <th>ClinGen SVI Calibrations</th>
      <td colspan="3">Abou Tayoun et al. 2018 (PVS1 LoF decision tree: NMD escape, last exon, &lt;10% CDS truncation, LOEUF &lt;0.60); ClinGen SVI 2020 (PM2_Supporting using non-bottlenecked reference population max AF across AFR, AMR, EAS, NFE, SAS, excluding founder populations ASJ, FIN, AMI, MID, REMAINING/OTH); ClinGen SVI PM3 (0.5 pts &rarr; PM3_Supporting for rare homozygous recessive variants); Pejaver et al. 2022 (graded PP3/BP4 REVEL and SpliceAI thresholds); Biesecker et al. 2018 (PP5/BP6 retired).</td>
    </tr>
    <tr>
      <th>Data &amp; Tool Versions</th>
      <td colspan="3">Reference: <code>GRCh38</code>; Pipeline: <code>{pipeline_version}</code>; Step image: <code>{image_tag}</code>; STRchive: <code>strchive_used={_esc(data.get("strchive_used"))}</code>.</td>
    </tr>
    <tr>
      <th>Workflow Limitations</th>
      <td colspan="3">Single-sample proband analysis without parental/pedigree genotypes: de novo (<code>PS2</code>, <code>PM6</code>) and co-segregation (<code>PP1</code>, <code>BS4</code>) criteria are <code>Not evaluable (phase_unknown)</code>. Curated functional assay databases are not wired (<code>PS3</code>, <code>BS3</code> <code>Not evaluable</code>). Repeat expansion sizing from short-read WGS has confidence intervals; clinical confirmation by PCR/Southern blot is recommended.</td>
    </tr>
    <tr>
      <th>Run Identifiers</th>
      <td colspan="3">Patient ID: <code>{patient_id}</code> (VCF sample: <code>{sample}</code>) | Schema version: <code>{_esc(data.get("schema_version", 1))}</code> | Pipeline: <code>{pipeline_version}</code> | Image: <code>{image_tag}</code> | LLM backend/model/status: <code>{llm_backend}</code> / <code>{llm_model}</code> / <code>{llm_status}</code>{(" — " + llm_error) if llm_error else ""}</td>
    </tr>
  </table>
  {appendix_c_html}

</body>
</html>
"""


def _import_weasyprint():
    try:
        import weasyprint  # type: ignore[import-untyped]
        return weasyprint
    except ModuleNotFoundError:
        root = Path(__file__).resolve().parent.parent
        for sp in root.glob(".venv/lib/python*/site-packages"):
            sp_str = str(sp)
            if sp_str not in sys.path:
                sys.path.append(sp_str)
        import weasyprint  # type: ignore[import-untyped]
        return weasyprint


def write_pdf(path: Path | str, data: dict) -> None:
    """Render ``data`` to a clinical PDF at ``path`` using WeasyPrint."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    html_str = render_pdf_html(data)
    weasyprint = _import_weasyprint()
    weasyprint.HTML(string=html_str).write_pdf(str(p))


def safe_write_pdf(path: Path | str, data: dict) -> bool:
    """Write ``path`` via WeasyPrint; on any error write ``<path>.error.txt`` and return False."""
    p = Path(path)
    try:
        write_pdf(p, data)
        return True
    except Exception as exc:
        p.parent.mkdir(parents=True, exist_ok=True)
        err_path = Path(f"{p}.error.txt")
        err_path.write_text(f"PDF generation failed for {p.name}: {type(exc).__name__}: {exc}\n")
        return False
