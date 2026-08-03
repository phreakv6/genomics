#!/usr/bin/env bash
# Download everything module 04 needs. Run this once, before run.sh.
#
#   ~45MB download, ~460MB on disk once indexed, ~3 minutes.
#   Safe to re-run: each step is skipped if its output already exists.
#
# Run: bash 04-pipeline/fetch_data.sh

set -euo pipefail

source "$(dirname "$0")/../check_env.sh"
require_tools bwa samtools curl

cd "$(dirname "$0")"
mkdir -p data && cd data

echo "== 1. chromosome 11 reference (GRCh38/hg38) ========================"
# The whole chromosome, not a slice, so that coordinates in our VCF are real
# genomic coordinates and no offset arithmetic is needed anywhere.
if [ ! -f chr11.fa ]; then
  echo "   downloading (~41MB compressed, 131MB expanded)..."
  curl -# -o chr11.fa.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr11.fa.gz
  gunzip -kf chr11.fa.gz
else
  echo "   chr11.fa already present, skipping"
fi
[ -f chr11.fa.fai ] || samtools faidx chr11.fa

echo
echo "== 2. bwa index (~90 seconds, ~400MB) =============================="
# Builds the Burrows-Wheeler/FM index from module 02. Done once per reference.
# For the full human genome this is ~1 hour and ~5GB -- one chromosome keeps
# it to something you can sit through.
if [ ! -f chr11.fa.bwt ]; then
  time bwa index chr11.fa
else
  echo "   index already present, skipping"
fi

echo
echo "== 3. real HG002 reads, chr11:5,200,000-5,300,000 =================="
# HG002 / NA24385: the GIAB Ashkenazi son, sequenced at 300x on Illumina.
#
# The BAM on NCBI's server is ~250GB. We do not download it. samtools reads
# the remote .bai index over HTTPS, works out which compressed blocks cover
# our region, and range-requests only those -- about 19MB. Remote indexed
# access is one of the most useful and least-known tricks in the field.
HG002=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam

if [ ! -f region.300x.bam ]; then
  echo "   pulling region out of a ~250GB remote BAM (downloads ~19MB)..."
  samtools view -b -o region.300x.bam "$HG002" chr11:5200000-5300000
  samtools index region.300x.bam
else
  echo "   region.300x.bam already present, skipping"
fi
echo "   reads: $(samtools view -c region.300x.bam)"

echo
echo "Done. Now run:  bash 04-pipeline/run.sh"
