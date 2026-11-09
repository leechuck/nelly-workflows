#!/usr/bin/env bash
# Stage the ClinVar GRCh38 VCF (+ index) for the ACMG agent from NCBI onto ws
# and print (NOT run) the arv-put command that uploads it to Keep as the
# `acmg_clinvar_vcf` workflow input. Refreshing ClinVar never requires an
# image rebuild; the gene table is built separately by build_gene_table.py.
#
# Requires: /mnt/data1/nelly-sim/acmg/clinvar.vcf.gz(.tbi) already present
# (downloaded by the supervisor; WP4 does not re-download on cluster nodes).
# Run only on ws (not on a laptop, not on the cluster, per COMMON rule 6).
set -euo pipefail

DIR=${ACMG_DATA_DIR:-/mnt/data1/nelly-sim/acmg}
VCF="$DIR/clinvar.vcf.gz"
TBI="$VCF.tbi"

if [ ! -f "$VCF" ] || [ ! -f "$TBI" ]; then
  echo "ClinVar VCF not found under $DIR; downloading GRCh38 ClinVar VCF..."
  mkdir -p "$DIR"
  cd "$DIR"
  wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
  wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_check_sum.md5
  md5sum -c clinvar_check_sum.md5 || echo "WARNING: md5 mismatch, verify manually"
fi

echo
echo "ClinVar files staged under $DIR:"
ls -l "$VCF" "$TBI"
echo
echo "Then upload to Keep (run by hand, NOT from this script):"
echo
echo "  arv-put --project <acmg-project-uuid> --progress"
echo "      --name 'ClinVar GRCh38 VCF (ACMG agent)'"
echo "      $VCF $TBI"
echo "  (one line above is the whole command without the indentation)"
echo "and copy the resulting collection PDH into the three main-vg*.json files as"
echo "  \"acmg_clinvar_vcf\": {\"class\": \"File\", \"location\": \"keep:<PDH>/clinvar.vcf.gz\","
echo "                        \"secondaryFiles\": [{\"class\": \"File\", \"location\": \"keep:<PDH>/clinvar.vcf.gz.tbi\"}]}."
