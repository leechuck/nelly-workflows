from pathlib import Path
import pytest
from acmg_agent.inputs import VepIndex, ClinvarIndex, load_gene_table, read_slivar_info
from acmg_agent.evidence import build_evidence
from acmg_agent.criteria import evaluate_all
from acmg_agent.classify import classify

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
        "path_score": "0.946", "pheno_score": "0.545", "combined": "0.9994",
        "disease": "Tay-Sachs disease (OMIM:272800)", "pheno_evidence": "",
    }
    ev = build_evidence(cand, variants_rows=[], vep=vep, clinvar=clinvar, gene_table=gt)
    assert ev.gene == "HEXA"
    assert ev.consequence == "splice_donor_variant"
    assert ev.mane_select is True
    assert ev.lof_mechanism is True
    assert ev.clinvar_sig == "Pathogenic"
    assert ev.clinvar_stars >= 2
    assert ev.max_af == pytest.approx(6.301e-6)
    assert ev.max_af_source == "gnomADe_NFE_AF"
    assert ev.max_af_any_pop == pytest.approx(0.003295)
    assert ev.max_af_any_pop_pop == "ASJ"
    assert ev.zygosity == "hom"
    assert ev.nmd_escape is False

    criteria = evaluate_all(ev)
    by_code = {c.code: c for c in criteria}
    assert (by_code["PM2"].status, by_code["PM2"].strength) == ("met", "supporting")
    assert "founder-population AF 0.0033 (ASJ)" in by_code["PM2"].rationale
    assert (by_code["PM3"].status, by_code["PM3"].strength) == ("met", "supporting")
    assert by_code["BS1"].status == "not_met"
    assert by_code["BA1"].status == "not_met"
    clf = classify(criteria)
    assert clf.points == 10
    assert clf.class_ == "Pathogenic"


def test_bottlenecked_populations_excluded_for_pm2_bs1_ba1(tmp_path):
    vcf_text = "\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr15,length=101991189>",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|SYMBOL|Feature|CANONICAL|MANE_SELECT|MAX_AF|gnomADe_AF|gnomADg_AF|gnomADe_ASJ_AF|gnomADe_FIN_AF|gnomADe_MID_AF|gnomADe_REMAINING_AF|gnomADg_AMI_AF|gnomADe_NFE_AF|gnomADg_EAS_AF">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
        "chr15\t72346234\t.\tC\tG\t65\tPASS\tCSQ=G|splice_donor_variant|HEXA|ENST00000268097.10|YES|NM_000520.6|0.06|0.002|0.001|0.06|0.04|0.02|0.03|0.05|1.0e-5|2.0e-5\tGT\t1/1",
    ]) + "\n"
    vcf_path = tmp_path / "founder.vcf"
    vcf_path.write_text(vcf_text)
    vep = VepIndex.from_vcf(vcf_path)
    clinvar = ClinvarIndex(FIX / "clinvar.vcf.gz")
    gt = load_gene_table(FIX / "gene_table.tsv")
    cand = {
        "gene": "HEXA", "moi": "AR", "chrom": "15", "pos": "72346234",
        "ref": "C", "alt": "G", "gt": "1/1", "effect": "splice_donor_variant",
        "hgvs": "HEXA:NM_000520.6:c.1421+1G>C:p.?", "rsid": "rs147324677",
        "clinvar": "PATHOGENIC", "max_freq": "",
        "path_score": "0.946", "pheno_score": "0.545", "combined": "0.9994",
        "disease": "Tay-Sachs disease (OMIM:272800)", "pheno_evidence": "",
    }
    ev = build_evidence(cand, [], vep, clinvar, gt)
    assert ev.max_af == pytest.approx(2.0e-5)
    assert ev.max_af_source == "gnomADg_EAS_AF"
    assert ev.max_af_any_pop == pytest.approx(0.06)
    assert ev.max_af_any_pop_pop == "ASJ"
    by_code = {c.code: c for c in evaluate_all(ev)}
    assert by_code["PM2"].status == "met"
    assert by_code["BS1"].status == "not_met"
    assert by_code["BA1"].status == "not_met"


def test_only_max_af_fallback_and_source_recorded():
    vep = VepIndex.from_vcf(FIX / "vep.vcf")
    clinvar = ClinvarIndex(FIX / "clinvar.vcf.gz")
    gt = load_gene_table(FIX / "gene_table.tsv")
    cand = {
        "gene": "ONLYMAXAF", "moi": "AD", "chrom": "22", "pos": "25000100",
        "ref": "C", "alt": "T", "gt": "0/1", "effect": "missense_variant",
        "hgvs": "", "rsid": "", "clinvar": "", "max_freq": "",
        "path_score": "", "pheno_score": "", "combined": "",
        "disease": "", "pheno_evidence": "",
    }
    ev = build_evidence(cand, [], vep, clinvar, gt)
    assert ev.max_af == pytest.approx(0.00042)
    assert ev.max_af_source == "MAX_AF"
    assert ev.max_af_any_pop == pytest.approx(0.00042)


def test_overall_gnomad_af_fallback_when_per_pop_absent():
    vep = VepIndex.from_vcf(FIX / "vep.vcf")
    clinvar = ClinvarIndex(FIX / "clinvar.vcf.gz")
    gt = load_gene_table(FIX / "gene_table.tsv")
    cand = {
        "gene": "AK3", "moi": "AD", "chrom": "9", "pos": "4719176",
        "ref": "C", "alt": "T", "gt": "0/1", "effect": "missense_variant",
        "hgvs": "", "rsid": "rs141473413", "clinvar": "", "max_freq": "",
        "path_score": "", "pheno_score": "", "combined": "",
        "disease": "", "pheno_evidence": "",
    }
    ev = build_evidence(cand, [], vep, clinvar, gt)
    assert ev.max_af == pytest.approx(0.000185)
    assert ev.max_af_source == "gnomADe_AF"


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
    assert ev.max_af is None
    assert ev.max_af_source is None


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
