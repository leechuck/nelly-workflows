"""Read-only evidence-pack tools and the one write tool for the ACMG agent.

Design §4.3 step 4: the model gets a bounded tool loop over read-only views
of the Evidence pack plus a single write tool ``submit_assessment``. Tools
never add data: every value returned already exists in the pack, which is
what makes the grounded-narrative check possible (agent.py).

Tool dispatch is a pure dict name -> callable(json_args_str) -> str, so the
agent loop stays tiny and the FakeBackend tests can exercise it directly.
"""
from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import asdict

from acmg_agent.model import Criterion, Evidence

READ_TOOL_NAMES = (
    "get_variant_evidence",
    "get_clinvar_same_variant",
    "get_clinvar_same_codon",
    "get_population_frequency",
    "get_gene_summary",
    "get_phenotype_match",
    "get_deterministic_criteria",
)
SUBMIT_TOOL = "submit_assessment"
TOOL_NAMES = (*READ_TOOL_NAMES, SUBMIT_TOOL)


def _read_tool_spec(name: str, description: str) -> dict:
    return {
        "type": "function",
        "function": {
            "name": name,
            "description": description,
            "parameters": {"type": "object", "properties": {}, "additionalProperties": False},
        },
    }


def read_tool_specs() -> list[dict]:
    """OpenAI-style function schemas for the read-only tools."""
    return [
        _read_tool_spec(
            "get_variant_evidence",
            "Variant identity and transcript annotation: chrom/pos/ref/alt, zygosity, "
            "MANE transcript, consequence, HGVS c./p., exon/intron and CDS/protein "
            "positions, NMD-escape flag, in-repeat flag, REVEL and SpliceAI scores.",
        ),
        _read_tool_spec(
            "get_clinvar_same_variant",
            "The ClinVar record for this exact variant (rsID/RCV id, interpretation, "
            "review stars). ClinVar labels are evidence for PS1/PM5 lookups only, "
            "never a classification shortcut.",
        ),
        _read_tool_spec(
            "get_clinvar_same_codon",
            "ClinVar P/LP lookups at this codon: max stars among same-amino-acid "
            "changes (PS1) and among different-amino-acid changes at the same codon "
            "(PM5).",
        ),
        _read_tool_spec(
            "get_population_frequency",
            "Population frequency data for the variant: max AF across sources, "
            "gnomAD popmax, homozygote count (nhomalt) and the gene-table BS1 "
            "threshold when available.",
        ),
        _read_tool_spec(
            "get_gene_summary",
            "Gene-level summary: MOI, gnomAD constraint (LOEUF, pLI), ClinVar "
            "missense/LoF P/LP composition, LoF-disease-mechanism flag, hotspot "
            "flag for this variant (not region details).",
        ),
        _read_tool_spec(
            "get_phenotype_match",
            "Phenotype fit: Exomiser human phenotype score, phenotype evidence "
            "string and disease name, plus any other variants in the same gene "
            "(their genotype, consequence and ClinVar interpretation).",
        ),
        _read_tool_spec(
            "get_deterministic_criteria",
            "The deterministic engine's 28-criterion assessment for this variant "
            "(code, status, strength, rationale). The agent adjusts these under "
            "the submission contract; it cannot introduce new codes.",
        ),
    ]


def all_tool_specs() -> list[dict]:
    """Read tools plus the single write tool, in the order sent to the model."""
    return [
        *read_tool_specs(),
        {
            "type": "function",
            "function": {
                "name": SUBMIT_TOOL,
                "description": (
                    "Submit the final assessment. criteria is a list of objects "
                    "{code, status, strength, reason} for every criterion you changed "
                    "or want to comment on (unchanged criteria may be omitted, they "
                    "keep the deterministic status). rationale is a short overall "
                    "justification grounded in the evidence pack."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "criteria": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "properties": {
                                    "code": {"type": "string"},
                                    "status": {"type": "string"},
                                    "strength": {"type": ["string", "null"]},
                                    "reason": {"type": "string"},
                                },
                                "required": ["code", "status"],
                            },
                        },
                        "rationale": {"type": "string"},
                    },
                    "required": ["criteria", "rationale"],
                },
            },
        },
    ]


class EvidencePack:
    """Immutable in-memory pack for one candidate variant.

    Everything the read tools may return is stored here, constructed once from
    the Evidence record and the deterministic criteria list.
    """

    def __init__(self, evidence: Evidence, det: list[Criterion]):
        self.evidence = evidence
        self.det = det

    # -- json renderings ---------------------------------------------------
    def as_json(self) -> str:
        return json.dumps(self.evidence_dict(), sort_keys=True)

    def evidence_dict(self) -> dict:
        d = asdict(self.evidence)
        d["zygosity"] = self.evidence.zygosity
        return d

    # -- read tools (name -> callable) -------------------------------------
    def dispatch(self) -> dict[str, Callable[[dict], str]]:
        return {
            "get_variant_evidence": self._tool_variant_evidence,
            "get_clinvar_same_variant": self._tool_clinvar_same_variant,
            "get_clinvar_same_codon": self._tool_clinvar_same_codon,
            "get_population_frequency": self._tool_population_frequency,
            "get_gene_summary": self._tool_gene_summary,
            "get_phenotype_match": self._tool_phenotype_match,
            "get_deterministic_criteria": self._tool_deterministic_criteria,
            SUBMIT_TOOL: self._tool_unsupported_submit,
        }

    def _tool_variant_evidence(self, args: dict) -> str:
        ev = self.evidence
        return json.dumps({
            "chrom": ev.chrom, "pos": ev.pos, "ref": ev.ref, "alt": ev.alt,
            "zygosity": ev.zygosity,
            "transcript": ev.transcript, "mane_select": ev.mane_select,
            "consequence": ev.consequence,
            "hgvsc": ev.hgvsc, "hgvsp": ev.hgvsp,
            "exon": ev.exon, "intron": ev.intron,
            "cds_position": ev.cds_position, "protein_position": ev.protein_position,
            "nmd_escape": ev.nmd_escape, "in_repeat_region": ev.in_repeat_region,
            "revel": ev.revel, "spliceai": ev.spliceai,
        }, sort_keys=True)

    def _tool_clinvar_same_variant(self, args: dict) -> str:
        ev = self.evidence
        return json.dumps({
            "clinvar_id": ev.clinvar_id, "clinvar_sig": ev.clinvar_sig,
            "clinvar_stars": ev.clinvar_stars,
        }, sort_keys=True)

    def _tool_clinvar_same_codon(self, args: dict) -> str:
        ev = self.evidence
        return json.dumps({
            "same_aa_ClinVar_PLP": ev.same_aa_ClinVar_PLP,   # max stars, PS1
            "diff_aa_ClinVar_PLP": ev.diff_aa_ClinVar_PLP,   # max stars, PM5
        }, sort_keys=True)

    def _tool_population_frequency(self, args: dict) -> str:
        ev = self.evidence
        return json.dumps({
            "max_af": ev.max_af, "popmax_af": ev.popmax_af,
            "nhomalt": ev.nhomalt, "bs1_threshold": ev.bs1_threshold,
        }, sort_keys=True)

    def _tool_gene_summary(self, args: dict) -> str:
        ev = self.evidence
        return json.dumps({
            "gene": ev.gene, "moi": ev.moi,
            "loeuf": ev.loeuf, "pli": ev.pli,
            "lof_mechanism": ev.lof_mechanism,
            "clinvar_missense_plp": ev.clinvar_missense_plp,
            "clinvar_missense_blb": ev.clinvar_missense_blb,
            "clinvar_lof_plp": ev.clinvar_lof_plp,
            "in_hotspot": ev.in_hotspot,
        }, sort_keys=True)

    def _tool_phenotype_match(self, args: dict) -> str:
        ev = self.evidence
        return json.dumps({
            "gene": ev.gene,
            "human_pheno_score": ev.human_pheno_score,
            "pheno_evidence": ev.pheno_evidence,
            "disease": ev.disease,
            "other_var_in_gene": ev.other_var_in_gene,
        }, sort_keys=True)

    def _tool_deterministic_criteria(self, args: dict) -> str:
        return json.dumps([
            {"code": c.code, "status": c.status, "strength": c.strength, "rationale": c.rationale}
            for c in self.det
        ], sort_keys=True)

    def _tool_unsupported_submit(self, args: dict) -> str:
        # submit_assessment must be routed through agent.py's contract
        # enforcement, never through the read-tool dispatch.
        return json.dumps({"error": "submit via the agent loop contract path"})


def build_prompt_context(pack: EvidencePack) -> str:
    """Compact context block for the first user message.

    The model may answer without tool calls when the pack is small; this
    serialisation is also the set of values the narrative grounding check
    recognises as 'present in the pack'.
    """
    return pack.as_json()
