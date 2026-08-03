#!/usr/bin/env bash
# The canonical FASTQ -> BAM -> VCF sequence, on a tiny known-truth dataset.
#
# Every human genomics pipeline in the world is this, plus more QC steps, more
# scale, and a heavier variant caller. Learn these six commands and the shape
# of the field is yours.
#
# Run: bash 03-formats/run.sh

set -euo pipefail

source "$(dirname "$0")/../check_env.sh"
require_tools bwa samtools bcftools

cd "$(dirname "$0")/data"

REF=HBB.hg38.fa
FQ=sample.hg38.fastq
SAMPLE=NA_sim

echo "== 1. index the reference =========================================="
# Builds the Burrows-Wheeler/FM index discussed in module 02. Done once per
# reference; for the whole human genome this takes ~1 hour and ~5GB.
bwa index "$REF" 2>&1 | sed 's/^/   /'
# .fai is a separate, simpler index letting tools seek to any coordinate.
samtools faidx "$REF"
ls -1 ${REF}* | sed 's/^/   /'

echo
echo "== 2. align reads -> SAM ==========================================="
# -R adds a read group (@RG). Variant callers refuse to run without one,
# because they need to know which sample each read came from.
bwa mem -R "@RG\tID:run1\tSM:${SAMPLE}\tPL:ILLUMINA" "$REF" "$FQ" \
    > aligned.sam 2> bwa.log
sed 's/^/   /' bwa.log | tail -3
echo "   -> aligned.sam ($(wc -l < aligned.sam | tr -d ' ') lines)"

echo
echo "== 3. sort and compress -> BAM ====================================="
# Coordinate-sorted BAM is the universal interchange format. Sorting is what
# makes indexing possible, and indexing is what makes random access possible.
samtools sort -o sorted.bam aligned.sam
samtools index sorted.bam
ls -lh aligned.sam sorted.bam sorted.bam.bai | awk '{print "   "$5"\t"$9}'

echo
echo "== 4. alignment statistics ========================================="
samtools flagstat sorted.bam | sed 's/^/   /'

echo
echo "== 5. call variants -> VCF ========================================="
# mpileup: for every position, gather the bases observed across all reads.
# call:    decide, per position, whether the pile is better explained by a
#          variant than by sequencing error -- a likelihood ratio using the
#          Phred quality scores from the FASTQ.
bcftools mpileup -f "$REF" sorted.bam -Ou 2>/dev/null \
  | bcftools call -mv -Ov -o calls.vcf 2>/dev/null
grep -vc '^#' calls.vcf | sed 's/^/   variant sites called: /'

echo
echo "== 6. the result ==================================================="
grep -v '^##' calls.vcf | sed 's/^/   /'

echo
echo "Ground truth planted by make_data.py: chr11:5227002 T>A, heterozygous."
echo "Our contig starts at chr11:5225464, so in contig coordinates that is"
echo "position $((5227002 - 5225464 + 1)). Check the POS column above."
