"""Clinical PDF report renderer (HTML + CSS -> PDF via WeasyPrint, design §11).

Main body presents the deterministic ACMG/AMP 2015 + ClinGen SVI evaluation:
  1. Header (patient/sample ID, sex, HPO terms + phenotype labels, GRCh38,
     pipeline/image versions, report date)
  2. Clinical Conclusion (leading finding, classification + points + combining
     rule, supported vs inconclusive decision, near-tie differential)
  3. ACMG/AMP 28-Criteria Evaluation Tables for the leading finding and every
     Pathogenic / Likely pathogenic candidate, grouped into Pathogenic (16) and
     Benign (12) criteria, showing decision status (Met + applied strength,
     Not met, Not evaluable + reason, Not applicable for retired PP5/BP6),
     evidence value/rationale, evidence source, point tally, and combining-rule
     result
  4. Ranked SNV / Indel Candidates Table
  5. Short Tandem Repeat Expansions (pathogenic/intermediate first, normal loci
     compact)
  Appendix A. Automated Reviewer Commentary (model, status, per-candidate
     adjustments + reasons, rationale, narrative — explicitly labelled as
     explanatory commentary, not the basis of the clinical decision)
  Appendix B. Methods, Data Versions, Limitations & Run Identifiers
"""
from __future__ import annotations

from datetime import datetime, timezone
import html
from pathlib import Path
import sys
from typing import Any

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


def _esc(val: Any, default: str = "—") -> str:
    if val is None or val == "":
        return default
    return html.escape(str(val))


def _criterion_source(code: str, ev: dict) -> str:
    code_u = (code or "").upper()
    if code_u == "PVS1":
        return "VEP 116 CSQ (MANE) + gene_table (LOEUF/LoF)"
    if code_u in ("PS1", "PM5"):
        return "ClinVar GRCh38 VCF + VEP 116 CSQ"
    if code_u in ("PM2", "BA1", "BS1"):
        src = ev.get("max_af_source")
        if src:
            return f"gnomAD v4 ({src})"
        return "VEP 116 gnomAD / slivar / Exomiser"
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


def _decision_label(status: str, strength: str | None, code: str) -> tuple[str, str, str]:
    """Return (decision_badge_html, applied_strength_str, css_row_class)."""
    st = (status or "").lower()
    is_benign = (code or "").upper().startswith("B")
    direction = "B" if is_benign else "P"
    if st == "met":
        pts = STRENGTH_POINTS.get((direction, strength or ""), 0)
        pts_str = f"+{pts}" if pts > 0 else str(pts)
        str_label = f"{strength or 'met'} ({pts_str} pt{'s' if abs(pts) != 1 else ''})" if strength != "stand_alone" else "stand_alone"
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


def _render_candidate_criteria_section(cand: dict, is_leading: bool = False) -> str:
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

    role_tag = "Leading SNV/Indel Candidate" if is_leading else "Pathogenic / Likely Pathogenic Candidate"
    conflict_note = (
        ' <span class="badge badge-warn">Pathogenic/Benign Evidence Conflict</span>'
        if conflict
        else ""
    )

    af_val = ev.get("max_af")
    af_str = f"{af_val:.2e} ({ev.get('max_af_source') or 'gnomAD'})" if isinstance(af_val, (int, float)) else "absent"
    if isinstance(ev.get("max_af_any_pop"), (int, float)) and (af_val is None or ev["max_af_any_pop"] > af_val):
        af_str += f" [founder pop {ev.get('max_af_any_pop_pop') or 'max'}: {ev['max_af_any_pop']:.4f}]"

    return f"""
    <div class="criteria-card">
      <h3>{role_tag}: {gene} — <code>{variant}</code></h3>
      <div class="cand-meta">
        <span><strong>HGVS:</strong> <code>{hgvs}</code></span>
        <span><strong>Zygosity / MOI:</strong> <code>{zyg}</code> / <code>{moi}</code></span>
        <span><strong>Disease:</strong> {disease}</span>
        <span><strong>Non-bottlenecked AF:</strong> <code>{_esc(af_str)}</code></span>
        <span><strong>Point Tally:</strong> <strong>{_esc(pts_str)} pts</strong></span>
        <span><strong>Richards Combining Rule:</strong> <strong>{rb_clf}</strong></span>
        <span><strong>Final Classification:</strong> <strong class="clf-label">{clf}</strong>{conflict_note}</span>
      </div>

      <h4>Pathogenic Criteria (PVS1, PS1–PS4, PM1–PM6, PP1–PP5)</h4>
      <table class="criteria-table">
        <thead>
          <tr>
            <th style="width: 7%;">Code</th>
            <th style="width: 21%;">Decision Status</th>
            <th style="width: 15%;">Applied Strength</th>
            <th style="width: 37%;">Evidence Value &amp; Evaluation Rationale</th>
            <th style="width: 20%;">Evidence Source</th>
          </tr>
        </thead>
        <tbody>
          {path_rows}
        </tbody>
      </table>

      <h4>Benign Criteria (BA1, BS1–BS4, BP1–BP7)</h4>
      <table class="criteria-table">
        <thead>
          <tr>
            <th style="width: 7%;">Code</th>
            <th style="width: 21%;">Decision Status</th>
            <th style="width: 15%;">Applied Strength</th>
            <th style="width: 37%;">Evidence Value &amp; Evaluation Rationale</th>
            <th style="width: 20%;">Evidence Source</th>
          </tr>
        </thead>
        <tbody>
          {benign_rows}
        </tbody>
      </table>
    </div>
    """


def render_pdf_html(data: dict) -> str:
    """Render the clinical report HTML+CSS document from assessment JSON data."""
    sample = _esc(data.get("sample") or "proband")
    sex = _esc(data.get("sex") or "unspecified")
    hpo_list = data.get("hpo_terms") or []
    hpo_str = ", ".join(str(t) for t in hpo_list) if hpo_list else "none provided"
    phenotype_labels = _esc(data.get("phenotype") or "—")
    report_date = _esc(data.get("report_date") or datetime.now(timezone.utc).strftime("%Y-%m-%d"))

    lead = data.get("leading_candidate")
    near_ties = data.get("near_tie_alternatives") or []
    cands = data.get("candidates") or []
    expansions = data.get("expansions") or []
    llm = data.get("llm") or {}
    llm_backend = _esc(llm.get("backend") or "none")
    llm_model = _esc(llm.get("model") or "none")
    llm_status = _esc(llm.get("status") or "disabled")

    # Determine whether conclusion is supported or inconclusive
    is_supported = False
    if isinstance(lead, dict):
        if lead.get("kind") == "expansion" and lead.get("status") == "pathogenic_range":
            is_supported = True
        elif lead.get("classification") in ("Pathogenic", "Likely pathogenic"):
            is_supported = True
    conclusion_status = "SUPPORTED" if is_supported else "INCONCLUSIVE"
    conclusion_badge_cls = "badge-met" if is_supported else "badge-warn"

    # Section 2: Conclusion block
    if isinstance(lead, dict) and lead.get("kind") == "expansion":
        lead_summary_html = f"""
        <table class="kv-table">
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

    diff_html = ""
    if near_ties:
        ties_str = ", ".join(
            f"<strong>{_esc(t.get('gene'))}</strong> (<code>{_esc(t.get('variant'))}</code>, {_esc(t.get('classification'))}, {_esc(t.get('points'))} pts)"
            for t in near_ties
        )
        diff_html = f"<p><strong>Near-tie differential candidates:</strong> {ties_str}</p>"
    else:
        diff_html = "<p><strong>Near-tie differential candidates:</strong> none within tie margin.</p>"

    # Section 3: 28-criteria tables for leading SNV (or top SNV) + all P/LP SNV candidates
    criteria_sections: list[str] = []
    shown_keys: set[tuple[str, str]] = set()
    lead_snv = None
    if isinstance(lead, dict) and lead.get("kind") == "snv":
        lead_snv = next(
            (c for c in cands if c.get("gene") == lead.get("gene") and c.get("variant") == lead.get("variant")),
            None,
        )
    if lead_snv is None and cands:
        lead_snv = cands[0]

    if lead_snv is not None:
        criteria_sections.append(_render_candidate_criteria_section(lead_snv, is_leading=True))
        shown_keys.add((str(lead_snv.get("gene")), str(lead_snv.get("variant"))))

    for c in cands:
        key = (str(c.get("gene")), str(c.get("variant")))
        if key in shown_keys:
            continue
        if c.get("classification") in ("Pathogenic", "Likely pathogenic"):
            criteria_sections.append(_render_candidate_criteria_section(c, is_leading=False))
            shown_keys.add(key)

    criteria_html = (
        "\n".join(criteria_sections)
        if criteria_sections
        else "<p><em>No SNV/indel candidates were evaluated against the 28 ACMG/AMP criteria.</em></p>"
    )

    # Section 4: Ranked candidates table
    if cands:
        cand_rows: list[str] = []
        for idx, c in enumerate(cands, 1):
            met_codes = [
                f"{cr['code']} ({cr['strength']})" if cr.get("strength") else cr["code"]
                for cr in (c.get("criteria") or [])
                if isinstance(cr, dict) and cr.get("status") == "met"
            ]
            cand_rows.append(
                f"<tr>"
                f"<td>{idx}</td>"
                f"<td><strong>{_esc(c.get('gene'))}</strong></td>"
                f"<td><code>{_esc(c.get('variant'))}</code></td>"
                f"<td><code>{_esc(c.get('hgvs'))}</code></td>"
                f"<td><code>{_esc(c.get('zygosity'))}</code></td>"
                f"<td><code>{_esc(c.get('moi'))}</code></td>"
                f"<td><strong>{_esc(c.get('classification'))}</strong></td>"
                f"<td>{_esc(c.get('points'))}</td>"
                f"<td>{_esc(', '.join(met_codes) if met_codes else '—')}</td>"
                f"<td>{_esc(c.get('disease'))}</td>"
                f"</tr>"
            )
        ranked_table_html = f"""
        <table class="data-table">
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
    else:
        ranked_table_html = "<p><em>No SNV/indel candidates evaluated.</em></p>"

    # Section 5: Repeat expansions (pathogenic/intermediate first, normal compact)
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

    # Appendix A: Automated reviewer commentary
    reviewed_genes = ", ".join(llm.get("reviewed") or []) or "none"
    unreviewed_genes = ", ".join(llm.get("unreviewed") or []) or "none"
    appendix_a_cards: list[str] = []

    if isinstance(lead, dict) and lead.get("kind") == "expansion":
        exp_st = _esc(lead.get("agent_status") or "not_reviewed")
        appendix_a_cards.append(f"""
        <div class="llm-card">
          <h4>{_esc(lead.get("gene"))} (Repeat Expansion <code>{_esc(lead.get("genotypes"))}</code>) — Status: <code>{exp_st}</code></h4>
          <ul>
            <li><strong>Review status:</strong> <code>{exp_st}</code></li>
            <li><strong>Adjusted ACMG criteria:</strong> none (repeat expansions are evaluated against locus thresholds, not SNV/indel ACMG criteria)</li>
            <li><strong>Reviewer rationale:</strong> {_esc(lead.get("agent_rationale"))}</li>
            <li><strong>Reviewer narrative:</strong> {_esc(lead.get("agent_narrative"))}</li>
          </ul>
        </div>
        """)

    for c in cands:
        cst = _esc(c.get("agent_status") or "not_reviewed")
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
        else:
            adj_str = "none"
        appendix_a_cards.append(f"""
        <div class="llm-card">
          <h4>{_esc(c.get("gene"))} (<code>{_esc(c.get("variant"))}</code>) — Status: <code>{cst}</code></h4>
          <ul>
            <li><strong>Review status:</strong> <code>{cst}</code></li>
            <li><strong>Adjusted criteria:</strong> {adj_str}</li>
            <li><strong>Reviewer rationale:</strong> {_esc(c.get("agent_rationale"))}</li>
            <li><strong>Reviewer narrative:</strong> {_esc(c.get("agent_narrative"))}</li>
          </ul>
        </div>
        """)

    appendix_a_body = (
        "\n".join(appendix_a_cards)
        if appendix_a_cards
        else "<p><em>No candidates evaluated.</em></p>"
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>ACMG/AMP Clinical Variant Report — {sample}</title>
  <style>
    @page {{
      size: A4;
      margin: 15mm 12mm 16mm 12mm;
      @bottom-left {{
        content: "ACMG/AMP Clinical Variant Interpretation Report — Sample: {sample} (GRCh38)";
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
    body {{
      font-family: "DejaVu Sans", sans-serif;
      font-size: 8.5pt;
      line-height: 1.35;
      color: #1a1a1a;
      margin: 0;
    }}
    code {{
      font-family: "DejaVu Sans Mono", monospace;
      font-size: 8pt;
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
    }}
    h3 {{
      font-size: 9.5pt;
      margin: 8px 0 4px 0;
      color: #1e293b;
    }}
    h4 {{
      font-size: 8.5pt;
      margin: 6px 0 3px 0;
      color: #334155;
    }}
    p {{
      margin: 4px 0;
    }}
    .kv-table, .data-table, .criteria-table {{
      width: 100%;
      border-collapse: collapse;
      margin: 4px 0 8px 0;
      font-size: 7.8pt;
    }}
    .kv-table th, .kv-table td,
    .data-table th, .data-table td,
    .criteria-table th, .criteria-table td {{
      border: 1px solid #cbd5e1;
      padding: 3px 5px;
      text-align: left;
      vertical-align: top;
    }}
    .kv-table th {{
      background-color: #f1f5f9;
      width: 18%;
      font-weight: bold;
    }}
    .kv-table td {{
      width: 32%;
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
    .cand-meta {{
      background-color: #f8fafc;
      border: 1px solid #cbd5e1;
      padding: 4px 6px;
      margin-bottom: 6px;
      font-size: 8pt;
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
    }}
    .llm-card ul {{
      margin: 3px 0;
      padding-left: 16px;
    }}
  </style>
</head>
<body>

  <h1>1. Clinical Variant Interpretation Report (ACMG/AMP 2015 + ClinGen SVI)</h1>
  <table class="kv-table">
    <tr>
      <th>Patient / Sample ID</th>
      <td><code>{sample}</code></td>
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
      <td colspan="3"><code>nelly-pangenome-vg</code> (<code>vg 1.78.0</code>, <code>DeepVariant 1.10.0</code>, <code>VEP 116</code>, <code>Exomiser 2512</code>, <code>nelly-acmg-agent:v4</code>)</td>
    </tr>
  </table>

  <h2>2. Clinical Conclusion</h2>
  {lead_summary_html}
  <p><strong>Summary:</strong> {_esc(data.get("narrative"))}</p>
  {diff_html}

  <h2>3. ACMG/AMP 28-Criteria Evaluation (Main Evidence Table)</h2>
  <p>Each of the 28 ACMG/AMP 2015 criteria is evaluated deterministically with ClinGen SVI calibrations. Every row states whether a decision was made (<strong>Met</strong> with applied strength, or <strong>Not met</strong>) or not made (<strong>Not evaluable</strong> with reason, or <strong>Not applicable</strong> for retired criteria PP5/BP6), together with the underlying evidence value, source, Tavtigian 2020 point tally, and Richards 2015 combining-rule classification.</p>
  {criteria_html}

  <h2>4. Ranked SNV / Indel Candidates</h2>
  {ranked_table_html}

  <h2>5. Short Tandem Repeat Expansions (ExpansionHunter + STRchive)</h2>
  {expansions_html}

  <div class="appendix-break"></div>
  <h2>Appendix A. Automated Reviewer Commentary ({llm_model}, {llm_status})</h2>
  <div class="notice-box">
    <strong>Explanatory Commentary Notice:</strong> This appendix contains the automated reviewer commentary (backend: <code>{llm_backend}</code>, model: <code>{llm_model}</code>, run status: <code>{llm_status}</code>). It is provided strictly as an explanation and audit trail, <strong>not as the basis of the clinical classification decision</strong>. All classifications in the main report are computed deterministically from the 28 ACMG/AMP criteria under fail-closed evidence-pack validation rules.
  </div>
  <p><strong>Reviewed genes:</strong> <code>{_esc(reviewed_genes)}</code> &nbsp;|&nbsp; <strong>Unreviewed genes:</strong> <code>{_esc(unreviewed_genes)}</code></p>
  {appendix_a_body}

  <h2>Appendix B. Methods, Data Versions, Limitations &amp; Run Identifiers</h2>
  <table class="kv-table">
    <tr>
      <th>ACMG/AMP Framework</th>
      <td colspan="3">Richards et al. 2015 (<em>Genet Med</em> 17:405) 28-criteria combining rules + Tavtigian et al. 2020 Bayesian point scale (PVS=8, PS=4, PM=2, PP=1, BS=&minus;4, BP=&minus;1; Pathogenic &ge;10 pts, Likely pathogenic 6–9 pts, VUS 0–5 pts, Likely benign &minus;1 to &minus;5 pts, Benign &le;&minus;6 pts or BA1).</td>
    </tr>
    <tr>
      <th>ClinGen SVI Calibrations</th>
      <td colspan="3">Abou Tayoun et al. 2018 (PVS1 LoF decision tree: NMD escape, last exon, &lt;10% CDS truncation, LOEUF &lt;0.60); ClinGen SVI 2020 (PM2_Supporting using non-bottlenecked reference population max AF across AFR, AMR, EAS, NFE, SAS, excluding founder populations ASJ, FIN, AMI, MID, REMAINING/OTH); ClinGen SVI PM3 (0.5 pts &rarr; PM3_Supporting for rare homozygous recessive variants); Pejaver et al. 2022 (graded PP3/BP4 REVEL and SpliceAI thresholds); Biesecker et al. 2018 (PP5/BP6 retired).</td>
    </tr>
    <tr>
      <th>Data &amp; Tool Versions</th>
      <td colspan="3">Reference: <code>GRCh38</code>; Alignment/Calling: <code>vg giraffe 1.78.0</code>, <code>DeepVariant 1.10.0</code>, <code>ExpansionHunter 5.0.0</code>, <code>bcftools 1.21</code>, <code>slivar 0.3.1</code>; Annotation: <code>Ensembl VEP 116</code> (gnomAD v4 exomes/genomes), <code>ClinVar GRCh38 VCF</code>, <code>Exomiser 15.0.0 (2512)</code>, <code>STRchive</code> (<code>strchive_used={_esc(data.get("strchive_used"))}</code>).</td>
    </tr>
    <tr>
      <th>Workflow Limitations</th>
      <td colspan="3">Single-sample proband analysis without parental/pedigree genotypes: de novo (<code>PS2</code>, <code>PM6</code>) and co-segregation (<code>PP1</code>, <code>BS4</code>) criteria are <code>Not evaluable (phase_unknown)</code>. Curated functional assay databases are not wired (<code>PS3</code>, <code>BS3</code> <code>Not evaluable</code>). Repeat expansion sizing from short-read WGS has confidence intervals; clinical confirmation by PCR/Southern blot is recommended.</td>
    </tr>
    <tr>
      <th>Run Identifiers</th>
      <td colspan="3">Sample: <code>{sample}</code> | Schema version: <code>{_esc(data.get("schema_version", 1))}</code> | Image: <code>nelly-acmg-agent:v4</code> | LLM backend/model/status: <code>{llm_backend}</code> / <code>{llm_model}</code> / <code>{llm_status}</code></td>
    </tr>
  </table>

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
