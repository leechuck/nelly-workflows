"""Core data structures for the ACMG/AMP 2015 + ClinGen SVI engine."""
from __future__ import annotations

from dataclasses import dataclass, field

CRITERION_STATUSES = ("met", "not_met", "not_evaluable", "not_applicable")


@dataclass
class Evidence:
    """Normalized variant-, gene-, population-, and case-level evidence pack."""

    gene: str
    moi: str                      # "AD" | "AR" | "XL" | "XR" | "XD" | "MT" | "OTHER"
    chrom: str
    pos: int
    ref: str
    alt: str
    zygosity: str                 # "het" | "hom" | "hemi" | "unknown"
    # VEP on MANE Select (or canonical fallback)
    transcript: str | None
    mane_select: bool
    consequence: str
    hgvsc: str | None
    hgvsp: str | None
    exon: str | None              # e.g. "11/14"
    intron: str | None            # e.g. "12/13"
    cds_position: int | None
    protein_position: int | None
    nmd_escape: bool              # True when stop/fs in last exon or last 50bp of penultimate
    in_repeat_region: bool        # True when inframe indel overlaps simple repeat
    # Population frequency
    max_af: float | None          # max of gnomADe_AF, gnomADg_AF, slivar popmax, Exomiser MAX_FREQ
    popmax_af: float | None
    nhomalt: int | None
    bs1_threshold: float | None
    # In-silico
    revel: float | None
    spliceai: float | None
    # ClinVar
    clinvar_id: str | None
    clinvar_sig: str | None       # "Pathogenic" | "Likely_pathogenic" | "Benign" | "Likely_benign" | "Uncertain_significance" | "Conflicting" | None
    clinvar_stars: int            # 0..4
    same_aa_ClinVar_PLP: int      # max review stars of a different-nt P/LP at the same AA -> PS1
    diff_aa_ClinVar_PLP: int      # max review stars of a different-AA P/LP at the same codon -> PM5
    # Gene-level
    lof_mechanism: bool | None
    loeuf: float | None
    pli: float | None
    clinvar_missense_plp: int | None
    clinvar_missense_blb: int | None
    in_hotspot: bool | None
    # Case-level (from Exomiser / slivar)
    other_var_in_gene: list[dict] = field(default_factory=list)   # for PM3 / BP2 / BP5
    human_pheno_score: float | None = None
    pheno_evidence: str = ""
    disease: str = ""
    clinvar_lof_plp: int | None = None
    alphamissense: float | None = None
    protein_length: int | None = None
    max_af_source: str | None = None
    max_af_any_pop: float | None = None
    max_af_any_pop_source: str | None = None
    max_af_any_pop_pop: str | None = None


@dataclass
class Criterion:
    """Single ACMG/AMP criterion outcome.

    Attributes:
        code: Criterion identifier ("PVS1" .. "BP7").
        status: One of ("met", "not_met", "not_evaluable", "not_applicable").
            "not_applicable" is used for retired criteria (PP5, BP6) and contributes 0 points.
        strength: One of ("stand_alone", "very_strong", "strong", "moderate", "supporting") or None.
        rationale: Human-readable explanation grounded in Evidence fields.
    """

    code: str                     # "PVS1" .. "BP7"
    status: str                   # "met" | "not_met" | "not_evaluable" | "not_applicable"
    strength: str | None          # "stand_alone" | "very_strong" | "strong" | "moderate" | "supporting" | None
    rationale: str


@dataclass
class Classification:
    """Combined 5-tier ACMG classification from points (primary) and Richards 2015 rules."""

    class_: str                   # "Pathogenic" | "Likely pathogenic" | "VUS" | "Likely benign" | "Benign"
    points: int
    rule_based_class: str
    conflict: bool
