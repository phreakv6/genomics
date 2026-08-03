# 04 — The real pipeline, and scoring yourself against truth

Module 03 was a rehearsal: we planted a variant and found it. Here we run the same shape of
pipeline on **real sequencing reads from a real person**, and then do the thing that makes this
module worth more than the other four combined — **measure how wrong we were**.

## The sample: HG002

Genome in a Bottle (GIAB) is a NIST-led consortium that took a handful of people, sequenced them on
every platform that exists (Illumina at 300×, PacBio, Oxford Nanopore, 10X, Complete Genomics,
SOLiD, Ion Torrent), ran many callers, and integrated the results into a **benchmark VCF**: a
best-available statement of what variants that person actually has.

HG002 (also called NA24385) is the son in an Ashkenazi Jewish trio — mother HG004, father HG003.
Having the parents matters: a variant in the child that appears in neither parent is either a rare
de novo mutation or, far more often, an error. Trio consistency is one of the strongest truth
signals available in genomics.

Two files matter:

- `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` — the truth variants.
- `..._benchmark_noinconsistent.bed` — the **high-confidence regions**. This is the crucial one.

## Why the high-confidence BED is the intellectually important file

GIAB does not claim to know HG002's whole genome. It claims to know ~90% of it. The other ~10% —
segmental duplications, long tandem repeats, the MHC, centromeres — is where all the platforms
disagreed, so GIAB withheld judgement and excluded those regions.

This means: **a variant you call outside the high-confidence BED cannot be scored.** It is not a
false positive; it is unknown. If you count it as an error you punish yourself for GIAB's ignorance;
if you count it as correct you flatter yourself. The only honest move is to restrict both your calls
and the truth set to the confident regions before comparing.

This is the same discipline as backtest survivorship bias. The benchmark defines the universe, and
evaluating outside the universe produces a number that means nothing. Most published "our caller
achieves 99.5%" claims live or die on exactly this restriction, and the field learned it the hard
way.

## The metrics

With truth (T) and your calls (C), restricted to confident regions:

- **True positives** — in both. **False positives** — in yours only. **False negatives** — in truth only.
- **Precision** = TP / (TP + FP) — of what I called, how much was real?
- **Recall** (sensitivity) = TP / (TP + FN) — of what was there, how much did I find?
- **F1** = harmonic mean of the two — the single number people quote.

The trade-off is a dial, not a fact about your pipeline: loosening the quality threshold raises
recall and lowers precision. A caller has a whole precision–recall *curve*, and quoting one point
without the threshold is how people mislead with benchmarks. Report SNPs and indels separately —
indels are far harder, and a combined number hides that.

Modern short-read pipelines get roughly: SNPs ~99.5% precision and recall; indels ~99% precision,
~95% recall, and much worse inside repeats. Our unpolished pipeline should land visibly below that,
and understanding *why* is the lesson.

## Variant representation: why comparison is harder than a set intersection

A deletion in a homopolymer can be written many equivalent ways:

```
ref  AAAAA        pos 100  REF=AAAA  ALT=AAA
                  pos 101  REF=AAAA  ALT=AAA
                  pos 103  REF=AA    ALT=A
```

All describe the same biological event. A naive position-and-string comparison calls two of these a
false positive *and* a false negative simultaneously — double-punishing a correct call. Fixes:

1. **Normalisation** (`bcftools norm`) — left-align indels and split multi-allelic sites into a
   canonical form. Cheap, necessary, and what we do here.
2. **`hap.py` / `vcfeval`** — the proper tools, which compare *haplotypes*: do the two VCFs, when
   applied to the reference, produce the same sequence? That's the semantically correct question and
   it's what GA4GH standardised on. Heavier to install; worth knowing exists.

Our comparison uses `bcftools norm` + `bcftools isec`, and the results should be read knowing that
some of the "errors" are representation artefacts rather than real mistakes. That caveat is itself
part of learning the field.

## What a real pipeline has that ours doesn't

We run: align → sort → call. Production adds, roughly in order of how much it matters:

- **Duplicate marking** (`samtools markdup`) — PCR amplification during library prep creates
  identical copies of one original fragment. Counting them as independent evidence turns one error
  into an apparent high-confidence variant. This one genuinely matters; we do it.
- **Base quality score recalibration** (BQSR) — sequencers' self-reported qualities are
  systematically biased by machine cycle and sequence context. GATK re-derives them empirically.
- **Local realignment / haplotype assembly** — around indels, per-read alignment is unreliable;
  GATK's HaplotypeCaller and freebayes reassemble the local region and call variants from candidate
  haplotypes rather than from a column-by-column pileup. This is the single biggest accuracy gap
  between `bcftools call` (what we use — a pileup caller) and GATK.
- **Filtering** — VQSR, or hard filters on depth, strand bias, mapping quality.

And the current state of the art is a fourth thing entirely: **DeepVariant** turns the pileup into
an image and runs a CNN over it. It beat hand-engineered statistical callers, which is a genuinely
interesting result about where the field's remaining edge lives — the biology was fine, the
likelihood model was the bottleneck.

## What to read

- Buffalo, *Bioinformatics Data Skills* — Ch. 11.
- The GIAB papers (Zook et al.) on how the benchmark was constructed. Short, and the methodology is
  the interesting part — it's a consensus-under-disagreement problem, not a biology problem.

## Then

`run.sh` pulls real HG002 reads for the beta-globin locus straight out of a remote indexed BAM,
realigns them, calls variants, and `benchmark.py` scores the result against the GIAB truth set
restricted to high-confidence regions.
