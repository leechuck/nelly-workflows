"""Candidate selection, evidence-tier ranking, and ExpansionHunter repeat evaluation.

Ports the full ranking and prompt-grounding logic from genome-linter-llm/blurb.py
verbatim so all existing ranking tests in ../tests/test_ranking.py pass when pointed
at acmg_agent.candidates, and also provides select_candidates / rank_key /
choose_leading / read_eh_expansions / table formatters for the ACMG engine.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import sys
import textwrap

# =====================================================================================
# Verbatim port from genome-linter-llm/blurb.py
# =====================================================================================

EH_PATHO_THRESHOLDS = {
    "HTT": (36, "Huntington disease", "autosomal dominant"),
    "ATXN1": (39, "Spinocerebellar ataxia type 1 (SCA1)", "autosomal dominant"),
    "ATXN2": (33, "Spinocerebellar ataxia type 2 (SCA2)", "autosomal dominant"),
    "ATXN3": (60, "Spinocerebellar ataxia type 3 / Machado-Joseph (SCA3)", "autosomal dominant"),
    "ATXN7": (37, "Spinocerebellar ataxia type 7 (SCA7)", "autosomal dominant"),
    "CACNA1A": (20, "Spinocerebellar ataxia type 6 (SCA6)", "autosomal dominant"),
    "ATXN8OS": (80, "Spinocerebellar ataxia type 8 (SCA8)", "autosomal dominant"),
    "ATXN10": (800, "Spinocerebellar ataxia type 10 (SCA10)", "autosomal dominant"),
    "TBP": (49, "Spinocerebellar ataxia type 17 (SCA17)", "autosomal dominant"),
    "PPP2R2B": (51, "Spinocerebellar ataxia type 12 (SCA12)", "autosomal dominant"),
    "ATN1": (48, "Dentatorubral-pallidoluysian atrophy (DRPLA)", "autosomal dominant"),
    "AR": (38, "Spinal and bulbar muscular atrophy (Kennedy disease)", "X-linked recessive"),
    "FXN": (66, "Friedreich ataxia", "autosomal recessive"),
    "FMR1": (200, "Fragile X syndrome (full mutation)", "X-linked"),
    "DMPK": (50, "Myotonic dystrophy type 1 (DM1)", "autosomal dominant"),
    "CNBP": (75, "Myotonic dystrophy type 2 (DM2)", "autosomal dominant"),
    "C9ORF72": (30, "ALS / frontotemporal dementia (C9orf72)", "autosomal dominant"),
    "CSTB": (30, "Progressive myoclonic epilepsy type 1 (EPM1)", "autosomal recessive"),
    "JPH3": (41, "Huntington disease-like 2 (HDL2)", "autosomal dominant"),
}


from acmg_agent.repeats import classify_repeat


def read_eh(path, strchive_loci=None):
    """Parse an ExpansionHunter JSON; return loci whose larger allele reaches the
    pathogenic threshold, as grounded evidence the LLM can cite."""
    flags = []
    if not path or not os.path.exists(path):
        return flags
    try:
        with open(path) as f:
            doc = json.load(f)
    except (ValueError, OSError) as e:
        sys.stderr.write("genome-linter-llm: WARN could not read EH json %s: %s\n" % (path, e))
        return flags
    for locus in (doc.get("LocusResults") or {}).values():
        for v in (locus.get("Variants") or {}).values():
            vid = v.get("VariantId") or ""
            if strchive_loci is not None:
                call = classify_repeat(v, strchive_loci)
                if call.tier == "pathogenic":
                    pmin = (call.ranges.get("pathogenic") or (0, None))[0] or 0
                    flags.append({
                        "gene": call.gene, "disease": call.disease, "moi": call.moi,
                        "genotype": call.genotype, "max_allele": call.max_allele, "threshold": pmin,
                        "unit": v.get("RepeatUnit") or "", "region": v.get("ReferenceRegion") or "",
                        "ci": call.ci,
                    })
                continue
            thr = EH_PATHO_THRESHOLDS.get(vid)
            if not thr:
                continue
            gt = str(v.get("Genotype") or "")
            alleles = [int(x) for x in gt.split("/") if x.strip().lstrip("-").isdigit()]
            if not alleles:
                continue
            maxa = max(alleles)
            threshold, disease, moi = thr
            if maxa >= threshold:
                flags.append({
                    "gene": vid, "disease": disease, "moi": moi,
                    "genotype": gt, "max_allele": maxa, "threshold": threshold,
                    "unit": v.get("RepeatUnit") or "", "region": v.get("ReferenceRegion") or "",
                    "ci": v.get("GenotypeConfidenceInterval") or "",
                })
    flags.sort(key=lambda f: f["max_allele"] - f["threshold"], reverse=True)
    return flags


def read_genes(path, top_n=None):
    """Exomiser genes TSV -> rows sorted by combined score (top_n=None keeps all)."""
    rows = []
    with open(path) as fh:
        r = csv.DictReader((l for l in fh), delimiter="\t")
        r.fieldnames = [(f or "").lstrip("#").strip() for f in (r.fieldnames or [])]
        for row in r:
            rows.append({(k or "").lstrip("#").strip(): v for k, v in row.items()})

    def score(row):
        for k in ("EXOMISER_GENE_COMBINED_SCORE", "COMBINED_SCORE", "ExomiserGeneCombinedScore"):
            if k in row and row[k] not in (None, "", "."):
                try:
                    return float(row[k])
                except ValueError:
                    pass
        return 0.0

    rows.sort(key=score, reverse=True)
    return rows if top_n is None else rows[:top_n]


def pick(row, *names, default=""):
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n]
    return default


def variants_for(path, symbols):
    """All Exomiser variant rows for the given gene symbols, keyed by symbol."""
    by_gene = {s: [] for s in symbols}
    if not os.path.exists(path):
        return by_gene
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        r.fieldnames = [(f or "").lstrip("#").strip() for f in (r.fieldnames or [])]
        for row in r:
            row = {(k or "").lstrip("#").strip(): v for k, v in row.items()}
            g = pick(row, "GENE_SYMBOL", "GENE")
            if g in by_gene:
                by_gene[g].append(row)
    return by_gene


TIER_REPEAT_EXPANSION = 0
TIER_BIALLELIC_RECESSIVE = 1
TIER_MONOALLELIC_DOMINANT = 2
TIER_OTHER = 3
TIER_LABEL = {
    TIER_REPEAT_EXPANSION: "pathogenic-range repeat expansion",
    TIER_BIALLELIC_RECESSIVE: "biallelic P/LP in recessive-mode gene, phenotype fit",
    TIER_MONOALLELIC_DOMINANT: "monoallelic P/LP in dominant-mode gene, phenotype fit",
    TIER_OTHER: "Exomiser order",
}
PHENO_FIT_MIN = 0.5
RECESSIVE_MOI = {"AR", "XR"}
DOMINANT_MOI = {"AD", "XD", "MT"}
ACMG_STRENGTH = {"PATHOGENIC": 3, "LIKELY_PATHOGENIC": 2, "UNCERTAIN_SIGNIFICANCE": 1}
PLP_MIN_STRENGTH = 2


def _float(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def gene_moi(row):
    moi = str(pick(row, "MOI", "MODE_OF_INHERITANCE", default="")).strip().upper()
    if not moi:
        m = re.search(r"_([A-Z]{2})$", str(pick(row, "ID", default="")))
        moi = m.group(1) if m else ""
    return moi


def zygosity(gt):
    alleles = [a for a in re.split(r"[/|]", str(gt or "").strip()) if a not in ("", ".")]
    if not alleles:
        return "unknown"
    if all(a == "0" for a in alleles):
        return "reference"
    if len(alleles) == 1:
        return "hemizygous"
    return "homozygous" if len(set(alleles)) == 1 else "heterozygous"


def acmg_class(vrow):
    return str(pick(vrow, "EXOMISER_ACMG_CLASSIFICATION", "ACMG_CLASSIFICATION",
                    default="")).strip().upper()


def is_contributing(vrow):
    return str(pick(vrow, "CONTRIBUTING_VARIANT", default="")).strip().lower() in ("1", "true")


def classify_gene(row, variants):
    sym = pick(row, "GENE_SYMBOL", "GENE", default="?")
    moi = gene_moi(row)
    mine = [v for v in variants if gene_moi(v) in ("", moi)]
    contributing = [v for v in mine if is_contributing(v)]
    plp = [v for v in contributing if ACMG_STRENGTH.get(acmg_class(v), 0) >= PLP_MIN_STRENGTH]

    pheno = _float(pick(row, "EXOMISER_GENE_PHENO_SCORE", "PHENO_SCORE", default="0"))
    human_raw = pick(row, "HUMAN_PHENO_SCORE", default=None)
    human = _float(human_raw) if human_raw is not None else None
    fit = pheno >= PHENO_FIT_MIN and (human is None or human > 0)

    zyg, acmg, tier, basis = "", "", TIER_OTHER, ""
    if plp:
        hom = [v for v in plp if zygosity(pick(v, "GENOTYPE", default="")) in ("homozygous", "hemizygous")]
        het_sites = {(pick(v, "CONTIG", default=""), pick(v, "START", "POS", default=""))
                     for v in plp if zygosity(pick(v, "GENOTYPE", default="")) == "heterozygous"}
        if hom:
            best = max(hom, key=lambda v: ACMG_STRENGTH.get(acmg_class(v), 0))
            zyg, acmg = zygosity(pick(best, "GENOTYPE", default="")), acmg_class(best)
            biallelic = True
        elif len(het_sites) >= 2:
            zyg = "compound heterozygous"
            acmg = min((acmg_class(v) for v in plp
                        if zygosity(pick(v, "GENOTYPE", default="")) == "heterozygous"),
                       key=lambda c: ACMG_STRENGTH.get(c, 0))
            biallelic = True
        else:
            best = max(plp, key=lambda v: ACMG_STRENGTH.get(acmg_class(v), 0))
            zyg, acmg = zygosity(pick(best, "GENOTYPE", default="")), acmg_class(best)
            biallelic = False

        if moi in RECESSIVE_MOI:
            if biallelic and fit:
                tier, basis = TIER_BIALLELIC_RECESSIVE, (
                    "%s %s variant(s) in a recessive-mode (%s) gene with a phenotype match"
                    % (zyg, acmg, moi))
            elif biallelic:
                basis = ("%s %s in a recessive-mode (%s) gene but no phenotype fit "
                         "(Exomiser phenotype evidence below its uninformative default or "
                         "no human-disease term match)" % (zyg, acmg, moi))
            else:
                basis = ("single heterozygous %s variant in a recessive-mode (%s) gene: "
                         "carrier state, not a diagnosis" % (acmg, moi))
        elif moi in DOMINANT_MOI:
            if fit:
                tier, basis = TIER_MONOALLELIC_DOMINANT, (
                    "%s %s variant in a dominant-mode (%s) gene with a phenotype match"
                    % (zyg, acmg, moi))
            else:
                basis = ("%s %s in a dominant-mode (%s) gene but no phenotype fit"
                         % (zyg, acmg, moi))
        else:
            basis = "%s %s variant, mode of inheritance %s" % (zyg, acmg, moi or "unknown")
    else:
        if contributing:
            best = max(contributing, key=lambda v: ACMG_STRENGTH.get(acmg_class(v), 0))
            zyg, acmg = zygosity(pick(best, "GENOTYPE", default="")), acmg_class(best) or "NOT_AVAILABLE"
        basis = "no pathogenic/likely-pathogenic contributing variant (%s)" % (acmg or "unclassified")

    return {
        "kind": "gene", "gene": sym, "tier": tier, "basis": basis,
        "moi": moi, "zygosity": zyg, "acmg": acmg, "fit": fit,
        "combined": _float(pick(row, "EXOMISER_GENE_COMBINED_SCORE", "COMBINED_SCORE", default="0")),
        "pheno": pheno, "human_pheno": human,
        "row": row, "variants": mine,
    }


def expansion_candidate(flag):
    return {
        "kind": "repeat_expansion", "gene": flag["gene"], "tier": TIER_REPEAT_EXPANSION,
        "basis": ("pathogenic-range repeat expansion: larger allele %d >= full-penetrance "
                  "threshold %d (%s, %s); not scoreable by Exomiser"
                  % (flag["max_allele"], flag["threshold"], flag["disease"], flag["moi"])),
        "moi": flag["moi"], "zygosity": "", "acmg": "", "fit": None,
        "combined": 0.0, "pheno": 0.0, "human_pheno": None,
        "excess": flag["max_allele"] - flag["threshold"], "eh": flag,
    }


def rank_candidates(genes, vmap, eh_flags=None, top_n=None):
    cands = [expansion_candidate(f) for f in (eh_flags or [])]
    for i, g in enumerate(genes):
        c = classify_gene(g, vmap.get(pick(g, "GENE_SYMBOL", "GENE", default="?"), []))
        c["exomiser_rank"] = i
        cands.append(c)

    def key(c):
        if c["kind"] == "repeat_expansion":
            return (c["tier"], -c["excess"], 0.0, 0)
        if c["tier"] == TIER_OTHER:
            return (c["tier"], 0, 0.0, c["exomiser_rank"])
        return (c["tier"], -ACMG_STRENGTH.get(c["acmg"], 0), -c["combined"], c["exomiser_rank"])

    cands.sort(key=key)
    return cands if top_n is None else cands[:top_n]


EH_LEAD_MARKER = "LEAD WITH THE REPEAT EXPANSION"
SINGLE_LEADER_MARKER = "single leading candidate"
TIE_MARKER = "effectively TIED"
INCONCLUSIVE_MARKER = "the result is inconclusive"
SUPPORTED_MARKER = "do NOT call the result inconclusive"
TIE_MARGIN = 0.05
MAX_VARIANTS_SHOWN = 4


def _candidate_lines(i, c):
    lines = []
    if c["kind"] == "repeat_expansion":
        f = c["eh"]
        lines.append("  %d. %s  REPEAT EXPANSION (ExpansionHunter): %s repeat genotype %s units "
                     "(larger allele %d >= pathogenic threshold %d), unit=%s at %s%s -> %s (%s)" % (
                         i, f["gene"], f["unit"] or "repeat", f["genotype"], f["max_allele"],
                         f["threshold"], f["unit"] or "?", f["region"] or "?",
                         (", CI " + f["ci"]) if f["ci"] else "", f["disease"], f["moi"]))
        lines.append("       rank basis: %s" % c["basis"])
        return lines
    g = c["row"]
    comb = pick(g, "EXOMISER_GENE_COMBINED_SCORE", "COMBINED_SCORE", default="?")
    phen = pick(g, "EXOMISER_GENE_PHENO_SCORE", "PHENO_SCORE", default="?")
    var = pick(g, "EXOMISER_GENE_VARIANT_SCORE", "VARIANT_SCORE", default="?")
    hum = pick(g, "HUMAN_PHENO_EVIDENCE", default="")
    lines.append("  %d. %s  combined=%s pheno=%s variant=%s MOI=%s genotype=%s ACMG=%s" %
                 (i, c["gene"], comb, phen, var, c["moi"] or "?", c["zygosity"] or "?",
                  c["acmg"] or "n/a"))
    lines.append("       rank basis: %s" % c["basis"])
    if hum:
        lines.append("       human-phenotype evidence: %s" % hum[:240])
    shown = sorted(c["variants"], key=lambda v: not is_contributing(v))[:MAX_VARIANTS_SHOWN]
    for v in shown:
        eff = pick(v, "VARIANT_EFFECT", "FUNCTIONAL_CLASS")
        path = pick(v, "EXOMISER_VARIANT_SCORE", "VARIANT_SCORE")
        freq = pick(v, "MAX_FREQUENCY", "MAX_FREQ", "FREQUENCY", default="0")
        cs = pick(v, "CLINVAR_PRIMARY_INTERPRETATION", "CLINVAR", default="")
        lines.append("       variant: %s genotype=%s effect=%s hgvs=%s acmg=%s pathScore=%s "
                     "maxFreq=%s clinvar=%s%s" %
                     (pick(v, "ID", "VARIANT", default="?"), pick(v, "GENOTYPE", default="?"),
                      eff, pick(v, "HGVS", default=""), acmg_class(v) or "n/a", path, freq, cs,
                      " [contributing]" if is_contributing(v) else ""))
    return lines


def leader_instruction(ranked):
    if not ranked:
        return "- No candidates are listed; state that the result is inconclusive."
    lead = ranked[0]
    if lead["kind"] == "repeat_expansion":
        exp = [c for c in ranked if c["kind"] == "repeat_expansion"]
        names = ", ".join(c["gene"] for c in exp)
        plural = len(exp) > 1
        return (
            "- %s: candidate%s %s %s a PATHOGENIC-RANGE repeat expansion%s. Exomiser cannot "
            "score repeat expansions, so %s absent from Exomiser's own gene scores by "
            "construction, not for lack of evidence. If the disease fits the patient's "
            "phenotype terms it is very likely THE diagnosis and you must lead with it, "
            "citing the gene, the repeat genotype and threshold, and the disease, and state "
            "explicitly whether the phenotype fits. %sDo NOT name any SNV/indel gene as the "
            "leading candidate; cover the SNV/indel candidates afterwards as the differential."
            % (EH_LEAD_MARKER, "s" if plural else "", names, "are" if plural else "is",
               "s" if plural else "", "they are" if plural else "it is",
               ("Present the %d expansions as a differential in the order listed. " % len(exp))
               if plural else ""))
    tied = [c for c in ranked
            if c["tier"] == lead["tier"] and c["acmg"] == lead["acmg"]
            and lead["combined"] - c["combined"] <= TIE_MARGIN]
    if len(tied) > 1:
        return ("- Candidates %s are %s (same rank basis, Exomiser combined scores within "
                "%.2f): present them as a differential, in the order listed, WITHOUT declaring "
                "one of them the single answer."
                % (", ".join(c["gene"] for c in tied), TIE_MARKER, TIE_MARGIN))
    return ("- Candidate 1 (%s) leads the field on the stated rank basis; you may name it as "
            "the %s, then cover the remaining candidates as the differential."
            % (lead["gene"], SINGLE_LEADER_MARKER))


def conclusion_instruction(ranked):
    if any(c["tier"] < TIER_OTHER for c in ranked):
        return ("- At least one candidate is supported (see its rank basis), so %s."
                % SUPPORTED_MARKER)
    return ("- No candidate has a pathogenic-range repeat expansion or a pathogenic/"
            "likely-pathogenic variant with a phenotype fit: state that %s."
            % INCONCLUSIVE_MARKER)


def build_prompt(phenotype, genes, vmap, eh_flags=None, ranked=None):
    ranked = rank_candidates(genes, vmap, eh_flags) if ranked is None else ranked
    lines = []
    lines.append("Patient HPO/phenotype terms: %s" % (phenotype or "(none provided)"))
    lines.append("")
    lines.append("Ranked candidates (pipeline ranking: ExpansionHunter pathogenic-range repeat "
                 "expansions and Exomiser SNV/indel genes, ordered by the rank basis given for "
                 "each candidate -- repeat expansion > biallelic P/LP in a recessive-mode gene "
                 "with phenotype fit > monoallelic P/LP in a dominant-mode gene with phenotype "
                 "fit > remaining genes in Exomiser combined-score order):")
    for i, c in enumerate(ranked, 1):
        lines.extend(_candidate_lines(i, c))
    if not eh_flags:
        lines.append("")
        lines.append("No pathogenic-range repeat expansion was reported by ExpansionHunter.")
    evidence = "\n".join(lines)

    instruction = textwrap.dedent("""
        You are a clinical genomics assistant. Using ONLY the evidence below, write a
        concise interpretation (5-10 sentences).

        Grounding rules (strict):
        - Discuss the candidates in the order listed and do not reorder them: the order is
          the pipeline's ranking. When you explain why a candidate ranks where it does,
          cite its stated rank basis; do not invent reasons of your own.
        %s
        - For each candidate you discuss, cite its scores, its genotype and ACMG class, the
          contributing variant(s) with their predicted effect and coordinates, and the mode
          of inheritance. For a repeat expansion cite the repeat genotype and threshold.
        - Do NOT introduce genes, diseases, or facts not in the evidence, and do NOT use
          outside clinical knowledge to argue for or against a candidate.
        %s

        --- EVIDENCE ---
        %s
        --- END EVIDENCE ---

        Interpretation:""").strip() % (leader_instruction(ranked),
                                       conclusion_instruction(ranked), evidence)
    return instruction


def run_llm(prompt, model_path):
    from llama_cpp import Llama
    llm = Llama(model_path=model_path, n_ctx=8192, n_threads=os.cpu_count(), verbose=False)
    out = llm.create_chat_completion(
        messages=[{"role": "user", "content": prompt}],
        temperature=0.2, max_tokens=600)
    return out["choices"][0]["message"]["content"].strip()


def interpret(prompt, model_path, no_llm=False):
    if no_llm:
        return prompt, "LLM skipped (GL_NO_LLM/--no-llm); showing grounded evidence only."
    try:
        return run_llm(prompt, model_path), ""
    except Exception as e:
        sys.stderr.write("genome-linter-llm: WARN LLM unavailable (%s); "
                         "falling back to grounded evidence.\n" % e)
        return prompt, "LLM unavailable (%s); showing grounded evidence only." % e


def write_report(output, ranked, eh_flags, text, note=""):
    with open(output, "w") as fh:
        fh.write("# Genome-linter interpretation (Exomiser + ExpansionHunter grounded)\n\n")
        if note:
            fh.write("_%s_\n\n" % note)
        if eh_flags:
            fh.write("Pathogenic-range repeat expansions (ExpansionHunter): " + "; ".join(
                "%s %s units (>=%d) -> %s" % (f["gene"], f["genotype"], f["threshold"], f["disease"])
                for f in eh_flags) + "\n\n")
        if ranked:
            fh.write("Ranked candidates (repeat expansions and Exomiser genes, re-ranked by "
                     "genotype/MOI/ACMG; see rank basis in the evidence): " + ", ".join(
                         "%d.%s%s" % (i, c["gene"],
                                      " (repeat expansion)" if c["kind"] == "repeat_expansion" else "")
                         for i, c in enumerate(ranked, 1)) + "\n\n")
            exo = [c["gene"] for c in sorted((c for c in ranked if c["kind"] == "gene"),
                                             key=lambda c: c["exomiser_rank"])]
            if exo:
                fh.write("Exomiser combined-score order (SNV/indel only): " + ", ".join(
                    "%d.%s" % (i, g) for i, g in enumerate(exo, 1)) + "\n\n")
        else:
            fh.write("Candidates: none returned.\n\n")
        fh.write(text + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--variants", default="")
    ap.add_argument("--phenotype", default=os.environ.get("GL_PHENOTYPE", ""))
    ap.add_argument("--output", required=True)
    ap.add_argument("--model", default=os.environ.get("GL_MODEL_PATH", "/model/model.gguf"))
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--no-llm", dest="no_llm",
                    default=os.environ.get("GL_NO_LLM", "") not in ("", "0", "false"),
                    action="store_true")
    ap.add_argument("--eh-json", dest="eh_json", default="")
    a = ap.parse_args()

    genes = read_genes(a.genes)
    symbols = [pick(g, "GENE_SYMBOL", "GENE", default="?") for g in genes]
    vmap = variants_for(a.variants, symbols) if a.variants else {}
    eh_flags = read_eh(a.eh_json)
    ranked = rank_candidates(genes, vmap, eh_flags, top_n=a.top)
    if not genes and not eh_flags:
        write_report(a.output, ranked, eh_flags,
                     "No candidate genes were returned by Exomiser and no pathogenic-range "
                     "repeat expansion was detected. The result is inconclusive; review the "
                     "upstream VCF and phenotype terms.")
        return

    prompt = build_prompt(a.phenotype, genes, vmap, eh_flags, ranked=ranked)
    text, note = interpret(prompt, a.model, no_llm=a.no_llm)
    write_report(a.output, ranked, eh_flags, text, note)


# =====================================================================================
# ACMG Engine Candidate Helpers (select_candidates, rank_key, choose_leading, tables)
# =====================================================================================

EH_LOCI = {
    "AFF2":     (34,  200, "XL", "FRAXE intellectual disability (OMIM:309548)"),
    "AR":       (34,  38,  "XL", "Spinal and bulbar muscular atrophy (OMIM:313200)"),
    "ATN1":     (35,  48,  "AD", "Dentatorubral-pallidoluysian atrophy (OMIM:125370)"),
    "ATXN1":    (35,  39,  "AD", "Spinocerebellar ataxia 1 (OMIM:164400)"),
    "ATXN10":   (32,  800, "AD", "Spinocerebellar ataxia 10 (OMIM:603516)"),
    "ATXN2":    (31,  33,  "AD", "Spinocerebellar ataxia 2 (OMIM:183090)"),
    "ATXN3":    (44,  60,  "AD", "Spinocerebellar ataxia 3 / Machado-Joseph (OMIM:109150)"),
    "ATXN7":    (19,  37,  "AD", "Spinocerebellar ataxia 7 (OMIM:164500)"),
    "ATXN8OS":  (50,  80,  "AD", "Spinocerebellar ataxia 8 (OMIM:608768)"),
    "C9ORF72":  (24,  60,  "AD", "FTD/ALS1 (OMIM:105550)"),
    "CACNA1A":  (18,  20,  "AD", "Spinocerebellar ataxia 6 (OMIM:183086)"),
    "CNBP":     (30,  75,  "AD", "Myotonic dystrophy 2 (OMIM:602668)"),
    "CSTB":     (3,   30,  "AR", "Unverricht-Lundborg progressive myoclonus epilepsy (OMIM:254800)"),
    "DIP2B":    (23,  250, "AD", "FRA12A intellectual disability (OMIM:136630)"),
    "DMPK":     (34,  50,  "AD", "Myotonic dystrophy 1 (OMIM:160900)"),
    "FMR1":     (44,  200, "XL", "Fragile X syndrome / FXTAS (OMIM:300624)"),
    "FXN":      (33,  66,  "AR", "Friedreich ataxia (OMIM:229300)"),
    "GLS":      (30,  100, "AR", "Glutaminase deficiency (OMIM:618339)"),
    "HTT":      (26,  36,  "AD", "Huntington disease (OMIM:143100)"),
    "JPH3":     (28,  41,  "AD", "Huntington disease-like 2 (OMIM:606438)"),
    "NIPA1":    (8,   10,  "AD", "Hereditary spastic paraplegia 6 / ALS risk (OMIM:600363)"),
    "NOP56":    (14,  650, "AD", "Spinocerebellar ataxia 36 (OMIM:614153)"),
    "NOTCH2NLC":(40,  60,  "AD", "Neuronal intranuclear inclusion disease (OMIM:603472)"),
    "PABPN1":   (10,  11,  "AD", "Oculopharyngeal muscular dystrophy (OMIM:164300)"),
    "PHOX2B":   (20,  24,  "AD", "Congenital central hypoventilation syndrome (OMIM:209880)"),
    "PPP2R2B":  (32,  51,  "AD", "Spinocerebellar ataxia 12 (OMIM:604326)"),
    "RFC1":     (11,  400, "AR", "CANVAS (OMIM:614575)"),
    "SAMD12":   (30,  400, "AD", "Familial adult myoclonic epilepsy 1 (OMIM:601068)"),
    "SOX3":     (15,  22,  "XL", "X-linked hypopituitarism (OMIM:300123)"),
    "STARD7":   (30,  200, "AD", "Familial adult myoclonic epilepsy 2 (OMIM:607876)"),
    "TBP":      (42,  49,  "AD", "Spinocerebellar ataxia 17 (OMIM:607136)"),
    "TCF4":     (40,  50,  "AD", "Fuchs endothelial corneal dystrophy 3 (OMIM:613267)"),
}

LOF_EFFECTS = {
    "frameshift_variant", "stop_gained", "splice_acceptor_variant",
    "splice_donor_variant", "start_lost", "stop_lost",
    "frameshift_elongation", "frameshift_truncation",
}


def fval(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def read_tsv(path: str | Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_eh_expansions(path: str | Path | None, strchive_loci: dict | None = None) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    with open(path) as f:
        data = json.load(f)
    results = []
    for locus_id, locus in data.get("LocusResults", {}).items():
        for var_id, var in locus.get("Variants", {}).items():
            gt = var.get("Genotype", "")
            ci = var.get("GenotypeConfidenceInterval", "")
            if not gt:
                continue
            alleles = [int(a) for a in gt.split("/") if a.isdigit()]
            if not alleles:
                continue
            if strchive_loci is not None:
                call = classify_repeat(var, strchive_loci)
                if call.tier == "unknown_locus":
                    continue
                status_map = {
                    "pathogenic": "pathogenic_range",
                    "intermediate": "intermediate_or_carrier",
                    "normal": "normal",
                }
                norm_max = (call.ranges.get("benign") or (None, 0))[1] or 0
                path_min = (call.ranges.get("pathogenic") or (0, None))[0] or 0
                results.append({
                    "locus_id": locus_id,
                    "gene": call.gene,
                    "genotypes": gt,
                    "ci": ci,
                    "allele1": alleles[0],
                    "allele2": alleles[-1],
                    "max_allele": call.max_allele,
                    "normal_max": norm_max,
                    "pathogenic_min": path_min,
                    "moi": call.moi,
                    "disease": call.disease,
                    "status": status_map.get(call.tier, "normal"),
                    "tier": call.tier,
                    "ranges": call.ranges,
                })
                continue
            if var_id != locus_id and var_id not in EH_LOCI:
                continue
            max_al = max(alleles)
            min_al = min(alleles)
            norm_max, path_min, moi, disease = EH_LOCI.get(
                locus_id, (40, 60, "?", "Repeat expansion locus")
            )
            eff_al = min_al if (moi == "AR" and len(alleles) == 2) else max_al
            if eff_al >= path_min:
                status = "pathogenic_range"
                tier_str = "pathogenic"
            elif max_al > norm_max:
                status = "intermediate_or_carrier"
                tier_str = "intermediate"
            else:
                status = "normal"
                tier_str = "normal"
            results.append({
                "locus_id": locus_id,
                "gene": locus_id,
                "genotypes": gt,
                "ci": ci,
                "allele1": alleles[0],
                "allele2": alleles[-1],
                "max_allele": max_al,
                "normal_max": norm_max,
                "pathogenic_min": path_min,
                "moi": moi,
                "disease": disease,
                "status": status,
                "tier": tier_str,
                "ranges": {
                    "benign": (None, norm_max),
                    "intermediate": (norm_max + 1, path_min - 1) if path_min > norm_max + 1 else (None, None),
                    "pathogenic": (path_min, None),
                },
            })
    order = {"pathogenic_range": 0, "intermediate_or_carrier": 1, "normal": 2}
    results.sort(key=lambda r: (order[r["status"]], -r["max_allele"], r["gene"]))
    return results


def select_candidates(genes: list[dict], variants: list[dict], top: int = 15) -> list[dict]:
    contributing = [v for v in variants if v.get("CONTRIBUTING_VARIANT") == "1"]
    by_gene_moi: dict[tuple[str, str], list[dict]] = {}
    for v in contributing:
        key = (v.get("GENE_SYMBOL", ""), v.get("MOI", ""))
        by_gene_moi.setdefault(key, []).append(v)

    seen_gene: set[str] = set()
    rows: list[dict] = []
    for g in genes:
        gene = g.get("#GENE_SYMBOL") or g.get("GENE_SYMBOL", "")
        moi = g.get("MOI", "")
        if not gene or gene in seen_gene:
            continue
        seen_gene.add(gene)
        vs = by_gene_moi.get((gene, moi), [])
        v0 = vs[0] if vs else {}
        pheno_ev = g.get("PHENO_EVIDENCE") or g.get("HUMAN_PHENO_EVIDENCE") or ""
        if g.get("MOI_DISEASE_MAP"):
            disease = g.get("MOI_DISEASE_MAP", "").split("$$")[0].split("->")[-1]
        elif v0.get("EXOMISER_ACMG_DISEASE_NAME"):
            did = v0.get("EXOMISER_ACMG_DISEASE_ID") or ""
            dname = v0.get("EXOMISER_ACMG_DISEASE_NAME") or ""
            disease = f"{dname} ({did})" if did else dname
        elif ":" in pheno_ev:
            disease = pheno_ev.split(":", 1)[0].strip()
        else:
            disease = ""

        path_score = (
            v0.get("EXOMISER_VARIANT_PATHOGENICITY_SCORE")
            or v0.get("EXOMISER_VARIANT_SCORE")
            or v0.get("MAX_PATH", "")
        )
        rows.append({
            "exomiser_rank": g.get("#RANK") or g.get("RANK", ""),
            "gene": gene,
            "moi": moi,
            "combined": g.get("EXOMISER_GENE_COMBINED_SCORE", ""),
            "pheno_score": g.get("EXOMISER_GENE_PHENO_SCORE", ""),
            "human_pheno_score": g.get("HUMAN_PHENO_SCORE", ""),
            "var_score": g.get("EXOMISER_GENE_VARIANT_SCORE", ""),
            "chrom": v0.get("CONTIG", ""),
            "pos": v0.get("START", ""),
            "ref": v0.get("REF", ""),
            "alt": v0.get("ALT", ""),
            "gt": v0.get("GENOTYPE", ""),
            "effect": v0.get("FUNCTIONAL_CLASS", ""),
            "hgvs": g.get("HGVS") or v0.get("HGVS", ""),
            "rsid": v0.get("RS_ID", ""),
            "clinvar": g.get("CLINVAR_PRIMARY_INTERPRETATION")
                       or v0.get("CLINVAR_PRIMARY_INTERPRETATION", ""),
            "max_freq": v0.get("MAX_FREQ", ""),
            "path_score": path_score,
            "disease": disease,
            "pheno_evidence": pheno_ev,
        })
        if len(rows) >= top:
            break
    return rows


def has_phenotype_fit(c: dict) -> bool:
    """True when Exomiser has a real human-disease phenotype match, not just a mouse/PPI hit."""
    pheno = fval(c.get("pheno_score"))
    ev = c.get("pheno_evidence", "")
    has_human_disease = bool(c.get("disease")) or "(OMIM:" in ev or "(ORPHA:" in ev
    return pheno >= 0.5 and has_human_disease


def rank_key(c: dict) -> tuple:
    clinvar = (c.get("clinvar") or "").upper()
    moi = c.get("moi", "")
    gt = c.get("gt", "")
    effect = (c.get("effect") or "").lower()
    af = fval(c.get("max_freq"), 0.0)
    path = fval(c.get("path_score"), 0.0)
    comb = fval(c.get("combined"), 0.0)
    pheno = fval(c.get("pheno_score"), 0.0)

    is_hom_or_comp = gt in ("1/1", "1|1") or "/" in gt and gt not in ("0/1", "0|1", "1|0")
    zyg_valid = (moi in ("AD", "XD", "XL", "MT")) or is_hom_or_comp
    rare = af <= 0.01
    is_path_clinvar = any(k in clinvar for k in ("PATHOGENIC", "LIKELY_PATHOGENIC")) \
        and "CONFLICTING" not in clinvar
    is_lof = any(e in effect for e in LOF_EFFECTS)
    pheno_fit = has_phenotype_fit(c)

    if not zyg_valid or not rare or (moi in ("AR", "XR") and is_hom_or_comp and not is_path_clinvar and not pheno_fit):
        tier = 5
    elif is_path_clinvar:
        tier = 0
    elif is_lof and path >= 0.8 and pheno_fit:
        tier = 1
    elif is_lof and path >= 0.8:
        tier = 2
    elif path >= 0.8 and pheno_fit:
        tier = 3
    else:
        tier = 4
    return (tier, -comb, -pheno, -path)


def choose_leading(candidates: list[dict], expansions: list[dict]) -> tuple[tuple | None, list[dict]]:
    path_exp = [e for e in expansions if e["status"] == "pathogenic_range"]
    ranked = sorted(candidates, key=rank_key) if candidates else []
    top_snv = ranked[0] if ranked else None
    top_tier = rank_key(top_snv)[0] if top_snv else 99

    if path_exp and top_tier > 0:
        return ("expansion", path_exp[0]), []

    if top_snv is None or top_tier >= 4:
        if path_exp:
            return ("expansion", path_exp[0]), []
        return None, []

    top_comb = fval(top_snv.get("combined"))
    alternatives = [
        c for c in ranked[1:5]
        if rank_key(c)[0] == top_tier and (top_comb - fval(c.get("combined"))) <= 0.05
    ]
    return ("snv", top_snv), alternatives


def format_candidates_table(rows: list[dict]) -> str:
    cols = [
        "#", "Gene", "MOI", "GT", "Variant", "Effect", "HGVS",
        "ClinVar", "MaxAF", "PathScore", "PhenoScore", "Combined", "Disease",
    ]
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for i, r in enumerate(rows, 1):
        var = f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}"
        if len(var) > 28:
            var = var[:25] + "..."
        out.append("| " + " | ".join([
            str(i), r["gene"], r["moi"], r["gt"], var,
            r["effect"], r["hgvs"][:40], r["clinvar"] or "-",
            r["max_freq"] or "0", r["path_score"], r["pheno_score"],
            r["combined"], r["disease"][:45],
        ]) + " |")
    return "\n".join(out)


def format_expansions_table(rows: list[dict]) -> str:
    if not rows:
        return "_No ExpansionHunter repeat catalog results provided._"
    cols = ["Gene", "Repeat alleles", "95% CI", "Normal ≤", "Pathogenic ≥", "MOI", "Status", "Disease"]
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in rows:
        out.append("| " + " | ".join([
            r["gene"], r["genotypes"], r["ci"],
            str(r["normal_max"]), str(r["pathogenic_min"]),
            r["moi"], r["status"], r["disease"],
        ]) + " |")
    return "\n".join(out)
