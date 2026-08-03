# 03 — Formats, coordinates, and the footguns

This is the module that has no intellectual glamour and causes 90% of real bugs. Buffalo's
*Bioinformatics Data Skills* is essentially a whole book about this layer, and it's the one that
was written for someone with your background.

## The format zoo

Each format exists at a specific stage of the pipeline. Learn them in pipeline order and they stop
looking arbitrary:

```
FASTA   the reference genome            "here is the sequence"
FASTQ   raw reads off the sequencer     "here is sequence + how confident we are"
SAM     aligned reads, text             "here is where each read goes on the reference"
BAM     aligned reads, binary+compressed same thing, 5x smaller, indexable
VCF     variants                        "here is where this sample differs from the reference"
BED/GFF annotation intervals            "this region is a gene / exon / promoter"
```

### FASTA — sequence

```
>chr11
ACGTACGTACGT...
```
A header line starting `>`, then sequence wrapped at 60–80 characters. The human reference is a
FASTA file, ~3GB. Wrapping is why naive line-by-line parsing breaks.

### FASTQ — reads with quality

Four lines per read, always:
```
@read_identifier
ACGTACGTACGT          <- the bases the machine called
+
IIIIIIIIIIII          <- how confident it was, one character per base
```

The quality line is the interesting part. Each character encodes a **Phred score**:

```
Q = -10 * log10(P(base is wrong))
```

So Q10 = 1 error in 10, Q20 = 1 in 100, Q30 = 1 in 1000. Q30 is the usual "good data" threshold.
The score is stored as ASCII with an offset of 33 (`chr(Q + 33)`), so `I` = ASCII 73 = Q40.

This matters more than it looks: every downstream statistical decision — is this a real variant or a
sequencing error? — is a likelihood computation that uses these numbers. Quality scores are the
reason variant calling is Bayesian rather than a counting exercise.

### SAM/BAM — alignments

One line per aligned read, 11 mandatory tab-separated columns. The ones that matter:

- `FLAG` — a **bitfield**. Bit 4 (`0x10`) = read is on the reverse strand, `0x4` = unmapped,
  `0x400` = PCR duplicate, and so on. A flag of `99` is not an ID, it's `1+2+32+64`. `samtools
  flags` decodes them.
- `RNAME`, `POS` — which chromosome, and the 1-based leftmost position.
- `MAPQ` — Phred-scaled confidence that the read is in the *right place*. MAPQ 0 means "this read
  maps equally well elsewhere" — usually a repeat. Filtering on MAPQ is one of the highest-leverage
  things you can do.
- `CIGAR` — a compact edit string: `100M` = 100 matched/mismatched bases; `50M2D48M` = 50 aligned,
  2-base deletion, 48 aligned; `10S90M` = 10 bases soft-clipped (present in the read, not aligned).
  Soft clips are a signal, not noise — piles of clipped reads at one position often mean a
  structural variant.

BAM is the same data BGZF-compressed and coordinate-sorted, with a `.bai` index so you can jump
straight to `chr11:5,227,002` without reading 50GB. Sort-then-index is a ritual you'll type
constantly. CRAM is the newer, smaller variant that stores only differences from the reference.

### VCF — variants

Header lines start `##`, a column header starts `#CHROM`, then one line per variant site:

```
#CHROM  POS      ID  REF  ALT  QUAL  FILTER  INFO         FORMAT   SAMPLE1
chr11   5227002  .   T    A    221   PASS    DP=35;AF=0.5 GT:DP:GQ 0/1:35:99
```

The genotype `GT` is the payload: `0/0` homozygous reference, `0/1` heterozygous, `1/1` homozygous
alternate. `|` instead of `/` means phased (you know which copy of the chromosome each allele is on).

The crucial conceptual point: **a VCF records differences from a reference, not your DNA.** It is a
diff, and it is meaningless without knowing which reference build it was made against.

### BED and GFF/GTF — intervals

BED is minimal: `chrom  start  end  [name  score  strand]`. GFF/GTF is the richer annotation format
used for gene models (which exons belong to which transcript). Both answer "what is at this
location", and `bedtools` is the interval-arithmetic swiss army knife over them.

## Footgun 1 — 0-based vs 1-based coordinates

This is the single most common source of silent off-by-one bugs in genomics, and it is unavoidable
because the field never standardised:

| Format | Coordinates |
|--------|-------------|
| BED | **0-based, half-open** `[start, end)` |
| GFF/GTF | 1-based, closed |
| SAM/BAM | 1-based (in the text; **0-based in the pysam API**) |
| VCF | 1-based |
| UCSC browser | 1-based display, 0-based files |

The first base of chromosome 1 is `1` in a VCF and `0` in a BED. A 10-base feature at the start of a
chromosome is `chr1 0 10` in BED and `chr1 1 10` in GFF. When you write a converter, that ±1 is
where you will lose an afternoon. `pysam` silently uses 0-based even though the file it's reading is
1-based — deliberately, to match Python slicing, and it catches everyone once.

## Footgun 2 — reference builds

`hg19`/`GRCh37` and `hg38`/`GRCh38` are different coordinate systems for the same genome. A position
in one is a *different place* in the other. Mixing them does not error — it silently produces
plausible-looking nonsense. Symptoms: your VCF's REF allele doesn't match the reference base at that
position, or annotations land in the wrong gene.

Compounding it, the naming conventions differ too: Ensembl calls it `11`, UCSC calls it `chr11`, and
tools will refuse to match them (or, worse, quietly find nothing). Half of real bioinformatics
debugging is chromosome-name and build reconciliation. `CrossMap`/`liftOver` convert between builds,
imperfectly — some regions have no equivalent.

Rule to adopt now: **record the build in every filename and every README.** `sample.hg38.bam`.

## Footgun 3 — Apple Silicon and bioconda

Relevant to this machine. Many bioconda packages have no `osx-arm64` build. The fix is to create the
environment under the Intel subdir and let Rosetta translate:

```
CONDA_SUBDIR=osx-64 conda create -n bio -c conda-forge -c bioconda samtools bcftools bwa ...
```

That's how this project's env was built. Without it, conda reports packages as simply
nonexistent, which is a confusing error for a package that obviously exists.

## The tools, and what each is for

| Tool | Job |
|------|-----|
| `samtools` | Everything BAM: view, sort, index, stats, depth, mpileup |
| `bcftools` | Everything VCF: call, filter, query, stats, isec (set operations) |
| `bedtools` | Interval arithmetic: intersect, merge, subtract, coverage |
| `bwa` / `minimap2` | Alignment (short reads / long reads) |
| `fastqc` | Read quality control |
| `GATK` | The heavyweight variant caller; more accurate, far more ceremony |
| `VEP` / `SnpEff` | Variant annotation — consequence, gene, protein change |
| `IGV` | Look at your BAM. Do this. Many bugs are visible and invisible to stats |

`samtools`/`bcftools`/`bedtools` compose over pipes exactly like Unix tools, because they were
designed by people who thought that way. Your shell instincts transfer directly.

## What to read

- Buffalo, *Bioinformatics Data Skills* — this module *is* that book. Ch. 10 (sequence formats),
  Ch. 11 (alignment/SAM), Ch. 9 (ranges/BED), and the reproducibility material.

## Then

`explore_formats.py` builds each format by hand from a real slice of the human genome and dissects
it field by field — including decoding a FLAG bitfield and a CIGAR string — then `coordinates.py`
demonstrates the 0-vs-1-based bug concretely, so you've seen it fail before it costs you a day.
