#!/usr/bin/env bash
# Real reads from a real person, realigned and called from scratch.
#
# Sample: HG002 / NA24385 (GIAB Ashkenazi son), Illumina 300x, GRCh38.
# Region: chr11:5,200,000-5,300,000 -- the beta-globin locus, including the
#         HBB gene from modules 01 and 03 and the upstream locus control region.
#
# Run: bash 04-pipeline/run.sh

set -euo pipefail

BIN=/opt/homebrew/Caskroom/miniconda/base/envs/bio/bin
export PATH="$BIN:$PATH"
cd "$(dirname "$0")/data"

REF=chr11.fa
REGION=chr11:5200000-5300000
EVAL_REGION=chr11:5205000-5295000     # trimmed: edges have partial coverage
SAMPLE=HG002
THREADS=8

echo "== 1. downsample 300x -> ~30x ======================================"
# 300x is a benchmark luxury, not what anyone sequences. 30x is the standard
# for a clinical/research whole genome, so that is the honest thing to test at.
# -s SEED.FRACTION : samtools reads the integer part as a seed.
samtools view -b -s 11.10 region.300x.bam > region.30x.bam
samtools index region.30x.bam
echo "   300x: $(samtools view -c region.300x.bam) reads"
echo "    30x: $(samtools view -c region.30x.bam) reads"

echo
echo "== 2. BAM -> FASTQ (back to raw, so we can realign honestly) ========"
# collate groups read pairs together by name; -n means name-based, which is
# what the fastq converter needs to emit R1/R2 in matched order.
samtools collate -@ $THREADS -u -O region.30x.bam tmp.collate \
  | samtools fastq -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null -n 2>/dev/null
echo "   R1: $(( $(wc -l < R1.fastq) / 4 )) reads"
echo "   R2: $(( $(wc -l < R2.fastq) / 4 )) reads"

echo
echo "== 3. align paired reads -> sorted BAM =============================="
# Paired-end is the real-world default: you sequence both ends of a ~350bp
# fragment. The known gap between mates is extra information -- it is how
# structural variants and repeats get resolved.
bwa mem -t $THREADS -R "@RG\tID:hiseq\tSM:${SAMPLE}\tPL:ILLUMINA" \
    $REF R1.fastq R2.fastq 2> bwa.log \
  | samtools fixmate -m -@ $THREADS - - \
  | samtools sort -@ $THREADS -o hg002.sorted.bam -
samtools index hg002.sorted.bam
tail -2 bwa.log | sed 's/^/   /'

echo
echo "== 4. mark PCR duplicates =========================================="
# Library prep amplifies fragments. Identical copies of one original molecule
# are not independent evidence -- counting them turns one early PCR error into
# an apparently high-confidence variant. This step genuinely moves the numbers.
samtools markdup -@ $THREADS -s hg002.sorted.bam hg002.markdup.bam 2> markdup.log
samtools index hg002.markdup.bam
grep -E "DUPLICATE TOTAL|EXCLUDED|READ" markdup.log | head -4 | sed 's/^/   /'

echo
echo "== 5. coverage check ==============================================="
samtools depth -r $EVAL_REGION hg002.markdup.bam \
  | awk '{s+=$3; n++} END {printf "   mean depth %.1fx over %d positions\n", s/n, n}'

echo
echo "== 6. call variants ================================================"
# -a AD,DP : ask for per-sample allele depths, which benchmark.py inspects.
# -m       : the multiallelic caller (better than the legacy -c consensus one)
# -v       : output variant sites only
bcftools mpileup -f $REF -r $REGION -a AD,DP -Ou hg002.markdup.bam 2>/dev/null \
  | bcftools call -mv -Oz -o calls.raw.vcf.gz 2>/dev/null
bcftools index -f calls.raw.vcf.gz
echo "   raw calls: $(bcftools view -H calls.raw.vcf.gz | wc -l | tr -d ' ')"

echo
echo "== 7. normalise ===================================================="
# Left-align indels and split multi-allelic records so that the same biological
# event is written the same way in both VCFs. Without this, comparison in
# step 9 double-counts representation differences as errors.
bcftools norm -f $REF -m -any -Oz -o calls.vcf.gz calls.raw.vcf.gz 2>/dev/null
bcftools index -f calls.vcf.gz
echo "   after norm: $(bcftools view -H calls.vcf.gz | wc -l | tr -d ' ')"

echo
echo "== 8. fetch the GIAB truth set ====================================="
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
# Pull only our region out of the remote indexed VCF -- no need for the 400MB file.
bcftools view -r $REGION "$GIAB/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" -Ou 2>/dev/null \
  | bcftools norm -f $REF -m -any -Oz -o truth.vcf.gz 2>/dev/null
bcftools index -f truth.vcf.gz
echo "   truth variants in region: $(bcftools view -H truth.vcf.gz | wc -l | tr -d ' ')"

if [ ! -f highconf.bed ]; then
  curl -s -o highconf.full.bed \
    "$GIAB/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
fi
# Restrict to our evaluation window.
printf "chr11\t5205000\t5295000\n" > eval.bed
bedtools intersect -a highconf.full.bed -b eval.bed > highconf.bed
awk '{s += $3 - $2} END {printf "   high-confidence bases in window: %d of 90000 (%.1f%%)\n", s, s/900}' highconf.bed

echo
echo "== 9. compare, restricted to high-confidence regions ==============="
bcftools view -R highconf.bed -Oz -o calls.hc.vcf.gz calls.vcf.gz 2>/dev/null
bcftools view -R highconf.bed -Oz -o truth.hc.vcf.gz truth.vcf.gz 2>/dev/null
bcftools index -f calls.hc.vcf.gz; bcftools index -f truth.hc.vcf.gz

rm -rf isec && bcftools isec -p isec truth.hc.vcf.gz calls.hc.vcf.gz 2>/dev/null
echo "   isec/0000.vcf = truth only  (false negatives): $(grep -vc '^#' isec/0000.vcf)"
echo "   isec/0001.vcf = calls only  (false positives): $(grep -vc '^#' isec/0001.vcf)"
echo "   isec/0002.vcf = both        (true positives):  $(grep -vc '^#' isec/0002.vcf)"

echo
echo "Now run:  python 04-pipeline/benchmark.py"
