#!/usr/bin/env bash
# Pull a slice of 1000 Genomes Phase 3, chromosome 22.
#
# Note what this does NOT do: download the 200MB chromosome file. bcftools
# reads the remote .tbi index, works out which compressed blocks cover the
# requested region, and range-requests only those. Remote indexed access is
# one of the most useful and least-known tricks in the field -- it works the
# same way against the GIAB BAM in module 04.
#
# Run: bash 05-popgen/fetch_data.sh

set -euo pipefail
source "$(dirname "$0")/../check_env.sh"
require_tools bcftools curl
cd "$(dirname "$0")"
mkdir -p data && cd data

KG=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
VCF=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# Who each sample is: population and superpopulation labels.
# These are used only for COLOURING the plot -- never fed to the PCA.
curl -s "$KG/integrated_call_samples_v3.20130502.ALL.panel" -o samples.panel
echo "samples in panel: $(( $(wc -l < samples.panel) - 1 ))"

# -r  region      5Mb of chr22 is plenty for global structure
# -m2 -M2 -v snps biallelic SNPs only (indels and multi-allelics complicate
#                 the 0/1/2 encoding for no gain here)
# -q 0.05:minor   minor allele frequency > 5%: common variants only
# annotate -x     drop INFO and all FORMAT fields except GT, purely for size
echo "fetching chr22:20-25Mb ..."
bcftools view -r 22:20000000-25000000 -m2 -M2 -v snps -q 0.05:minor \
    "$KG/$VCF" -Ou 2>/dev/null \
  | bcftools annotate -x INFO,^FORMAT/GT -Oz -o chr22.slice.vcf.gz 2>/dev/null
bcftools index -f chr22.slice.vcf.gz

echo "variants: $(bcftools view -H chr22.slice.vcf.gz | wc -l | tr -d ' ')"
echo "samples:  $(bcftools query -l chr22.slice.vcf.gz | wc -l | tr -d ' ')"
ls -lh chr22.slice.vcf.gz

echo
echo "Note this is GRCh37 and uses Ensembl-style chromosome names ('22', not"
echo "'chr22') -- unlike the GRCh38/'chr11' data in module 04. Mixing the two"
echo "silently produces zero overlaps. See 03-formats/theory.md, footgun #2."
