from pathlib import Path
import pytest
from acmg_agent.inputs import VepIndex, ClinvarIndex, load_gene_table, read_slivar_info
from acmg_agent.evidence import build_evidence

FIX = Path(__file__).parent / "fixtures" / "acmg"


def test_hexa_evidence_fields():
    vep = VepIndex.from_vcf(FIX / "vep.vcf")
    clinvar = ClinvarIndex(FIX / "clinvar.vcf.gz")
    gt = load_gene_table(FIX / "gene_table.tsv")
    cand = {
        "gene": "HEXA", "moi": "AR", "chrom": "15", "pos": "72346234",
        "ref": "C", "alt": "G", "gt": "1/1", "effect": "splice_donor_variant",
        "hgvs": "HEXA:NM_000520.6:c.1421+1G>C:p.?", "rsid": "rs147324677",
        "clinvar": "PATHOGENIC", "max_freq": "6.3013897E-4",
        "path_score": "0.946", "pheno_score": "0.9815", "combined": "0.9994",
        "disease": "Tay-Sachs disease (OMIM:272800)", "pheno_evidence": "",
    }
    ev = build_evidence(cand, variants_rows=[], vep=vep, clinvar=clinvar, gene_table=gt)
    assert ev.gene == "HEXA"
    assert ev.consequence == "splice_donor_variant"
    assert ev.mane_select is True
    assert ev.lof_mechanism is True
    assert ev.clinvar_sig == "Pathogenic"
    assert ev.clinvar_stars >= 2
    assert ev.max_af == pytest.approx(6.3013897e-6)
    assert ev.zygosity == "hom"
    assert ev.nmd_escape is False


def test_missing_annotations_default_none():
    vep = VepIndex.from_vcf(FIX / "vep.vcf")
    clinvar = ClinvarIndex(FIX / "clinvar.vcf.gz")
    gt = load_gene_table(FIX / "gene_table.tsv")
    cand = {
        "gene": "UNKNOWNGENE", "moi": "AD", "chrom": "1", "pos": "999999999",
        "ref": "A", "alt": "T", "gt": "0/1", "effect": "missense_variant",
        "hgvs": "", "rsid": "", "clinvar": "", "max_freq": "",
        "path_score": "", "pheno_score": "", "combined": "",
        "disease": "", "pheno_evidence": "",
    }
    ev = build_evidence(cand, [], vep, clinvar, gt)
    assert ev.revel is None
    assert ev.spliceai is None
    assert ev.clinvar_sig is None
    assert ev.lof_mechanism is None


def test_same_codon_lookup():
    clinvar = ClinvarIndex(FIX / "clinvar.vcf.gz")
    hit = clinvar.same_codon("15", 72346234, "C", "G", hgvsp="p.Arg170Trp")
    assert hit is None or isinstance(hit, tuple)
    cftr_same = clinvar.same_codon("7", 117559590, "C", "T", hgvsp="ENSP00000003084.6:p.Ser508Phe")
    assert cftr_same == (True, 3)
    cftr_diff = clinvar.same_codon("7", 117559590, "C", "T", hgvsp="ENSP00000003084.6:p.Ser508Cys")
    assert cftr_diff == (False, 3)


def test_slivar_info_merged():
    vep = VepIndex.from_vcf(FIX / "vep.vcf")
    clinvar = ClinvarIndex(FIX / "clinvar.vcf.gz")
    gt = load_gene_table(FIX / "gene_table.tsv")
    slivar = read_slivar_info(FIX / "slivar.vcf")
    cand = {
        "gene": "HEXA", "moi": "AR", "chrom": "15", "pos": "72346234",
        "ref": "C", "alt": "G", "gt": "1/1", "effect": "splice_donor_variant",
        "hgvs": "HEXA:NM_000520.6:c.1421+1G>C:p.?", "rsid": "rs147324677",
        "clinvar": "PATHOGENIC", "max_freq": "6.3013897E-4",
        "path_score": "0.946", "pheno_score": "0.9815", "combined": "0.9994",
        "disease": "Tay-Sachs disease (OMIM:272800)", "pheno_evidence": "",
    }
    ev = build_evidence(cand, [], vep, clinvar, gt, slivar_info=slivar)
    assert ev.nhomalt == 0
    assert ev.popmax_af == pytest.approx(6.3013897e-6)
