# Chapter 10 — Variants on Disk: VCF

The last stop in the pipeline. A **VCF** (Variant Call Format) file records
where a sample differs from the reference — and *only* where it differs.
This chapter dissects the format, watches a variant caller produce one from
the module's BAM, and closes with the two conceptual points that people
misread VCFs by: a VCF is a diff, and absence is not evidence of absence.

## Calling variants: what mpileup actually does

Module 03's `run.sh`, step 5:

```bash
# mpileup: for every position, gather the bases observed across all reads.
# call:    decide, per position, whether the pile is better explained by a
#          variant than by sequencing error -- a likelihood ratio using the
#          Phred quality scores from the FASTQ.
bcftools mpileup -f HBB.hg38.fa sorted.bam -Ou | bcftools call -mv -Ov -o calls.vcf
```

The **pileup** is the pivotal data structure of variant calling: walk the
sorted BAM position by position, and at each position stack up every base
from every read covering it. `dissect.py` renders the pileup at rs334 —
this is what the caller sees at the planted site:

```text
The pileup -- what the caller actually sees
--------------------------------------------------------------------------
  chr11:5,227,002  reference base T
  38 reads stacked here:
    A  ######################                    22  (alternate)
    T  ################                          16  (reference)

  Two bases, roughly 50/50, both at high quality. That is the
  signature of a heterozygous site. One base at 50/50 with LOW
  qualities would be noise; a 5% minority would be contamination
  or somatic mosaicism. The caller distinguishes these purely
  from the Phred scores in the FASTQ.
```

This is where Part I and Part III snap together. The sample carries one
reference `T` and one sickle `A` (heterozygous, Chapter 4), so reads drawn
from the two chromosome copies split ~50/50 — 16 vs 22 here, binomial
noise around half. The caller's question at every position is a model
comparison: *is this stack better explained by genotype T/T plus
sequencing error, T/A, or A/A?* The Phred scores (Chapter 8) are the error
model's parameters; with 22 high-quality `A`s, the probability that T/T
plus independent Q30+ errors produced this stack is astronomically small,
and T/A wins. That likelihood computation — priors, error rates, posterior
per genotype — is why variant calling is Bayesian rather than "count
bases, apply threshold".

## Anatomy of a VCF record

The output file has `##` header lines (declaring the reference, the
caller, every field's definition), a `#CHROM` column-header line, then one
line per variant site. The module's entire result — one line, which is the
correct answer:

```text
#CHROM                 POS   ID  REF  ALT  QUAL     FILTER  INFO             FORMAT     NA_sim
chr11:5225464-5229395  1539  .   T    A    222.283  .       DP=38;MQ=60;...  GT:PL:AD   0/1:255,0,240:16,22
```

Dissected, with the numbers explained:

```text
  POS     1539   (contig-local; real coordinate chr11:5,227,002)
  REF     T      <- what the reference genome says
  ALT     A      <- what this sample has
  QUAL    222.3  (Phred: P(no variant here) = 5.91e-23)
  INFO    DP=38  <- 38 reads covered this position
  SAMPLE  NA_sim
    GT    0/1   <- 0=REF allele, 1=ALT allele. heterozygous
    AD    16,22  <- 16 reads showed T, 22 showed A
```

Field by field:

- **CHROM/POS** — where, **1-based**. (Here POS is local to our little
  contig, whose name records that it starts at chr11:5,225,464 — so the
  real coordinate is 5225464 + 1539 − 1 = 5,227,002. Off-by-one arithmetic
  even in the sanity check; welcome to Chapter 11.)
- **REF/ALT** — the reference allele and the alternate the sample carries.
- **QUAL** — Phred-scaled (the scale's third appearance) confidence that
  *some* variant is here at all: 222 ⇒ P(site is actually reference) ≈
  10^(−22.2).
- **GT**, the genotype, is the payload: `0/0` homozygous reference, `0/1`
  heterozygous, `1/1` homozygous alternate — Chapter 4's zygosity
  vocabulary, now a data field. (A `|` separator instead of `/` means
  **phased**: you additionally know which parental chromosome copy carries
  which allele.)
- **AD**, allelic depth, is the pileup's summary: 16 reference reads, 22
  alternate — the evidence behind the `0/1`.

And the header's closing observation, which ties three chapters together:

```text
  This is rs334. Reference says T at chr11:5,227,002;
  this sample carries one T and one A. HBB is on the minus strand,
  so on the coding strand that is the A>T that module 01 turned into
  p.Glu7Val -- sickle cell. Same event, three coordinate conventions.
```

Genome plus strand `T>A`, coding strand `A>T`, protein `p.Glu7Val`: one
physical event, three correct names, and you now hold the conversion rules
between all of them.

## The two ways people misread a VCF

**A VCF is a diff, not a genome.** It records differences *from a
reference build*, and it is meaningless without knowing which one. The
same physical variant is `chr11:5227002` against hg38 and `chr11:5248232`
against hg19 (Chapter 11); a VCF with no declared build is not data. This
is also why REF is in the file at all — it's redundant with the reference,
and that redundancy is the integrity check: if the VCF's REF doesn't match
the reference you think it was called against, you've caught a build
mismatch before it cost you a week.

**Absence means "reference", not "no data" — except when it doesn't.** A
standard VCF contains only sites where a variant was called. A position
missing from the VCF might be confidently reference — or might have had no
coverage at all, and the file looks identical either way. Distinguishing
"ref" from "no evidence" requires going back to the BAM (is there depth
there?) or emitting gVCF-style reference blocks. Part IV's benchmark
analysis ends on exactly this misreading — it checks whether HG002's VCF
"contains" the sickle site, and the answer teaches the lesson.

## Verifying against ground truth

Module 03's punchline: `make_data.py` planted exactly one variant, and the
pipeline's VCF contains exactly one record — the right one, with the right
genotype, and zero false positives anywhere else in the 3,932bp region.
The final lines of `run.sh` do the coordinate arithmetic and invite you to
check the POS column yourself.

That closes the loop this Part opened: reference in (FASTA), reads in
(FASTQ), alignments through (SAM/BAM), one heterozygous `T>A` out (VCF),
truth recovered. Everything Part IV does is this loop again — with real
reads, a real truth set, ninety thousand positions, and honest scoring.

> **BED and GFF, the formats this figure keeps gesturing at.** The
> remaining families in Figure 8-1 are *interval* formats — they answer
> "what is at this location" rather than "what is the sequence/read/
> variant". BED is the minimal one (`chrom start end`, 0-based half-open)
> and appears throughout Part IV, where a BED file defines the regions
> where truth is known. GFF/GTF is the rich one — gene models: which exons
> belong to which transcript on which strand, Chapter 4's "annotation"
> made file. It gets Chapter 12 to itself.

## Run it

```bash
bash 03-formats/run.sh          # steps 5-6 are this chapter
python 03-formats/dissect.py    # VCF + pileup sections
```

## Further reading

- Buffalo, *Bioinformatics Data Skills* — the VCF material.
- The [VCF specification](https://samtools.github.io/hts-specs/) — the
  INFO/FORMAT vocabulary is large; the spec is where it's defined.
