# Glossary

Terms in the order you'd want them when one has scrolled out of memory —
alphabetical, each with a pointer to the chapter that introduces it
properly.

**Allele** — one of the alternative versions of the sequence at a
position. You carry two (one per parental chromosome copy). *(Ch. 4)*

**Alternative splicing** — one gene spliced multiple ways, producing
different proteins from the same locus; the reason variant consequences
are per-transcript. *(Ch. 3, 4)*

**Amino acid** — one of the 20 building blocks of proteins; each has a
one-letter (`E`) and three-letter (`Glu`) code and a chemical personality
(charged, greasy, small…). *(Ch. 2)*

**Annotation** — curated knowledge about where features (genes, exons,
transcripts) sit on the genome, stored in GFF/GTF files; not derivable
from the sequence itself. *(Ch. 3, 4, 12)*

**BAM** — binary, BGZF-compressed, indexable SAM; the standard container
for aligned reads. *(Ch. 9)*

**Base / base pair (bp)** — one letter of DNA (`A`, `C`, `G`, `T`); bp is
the unit of sequence length (kb = 10³, Mb = 10⁶, Gb = 10⁹). *(Ch. 1)*

**BED** — minimal interval format: `chrom start end`, **0-based
half-open**; used to define regions (e.g. GIAB's high-confidence
universe). *(Ch. 10, 11, 13)*

**BLOSUM62** — the standard amino-acid substitution matrix; log-odds
scores from observed substitutions — an empirical prior for protein
alignment. *(Ch. 6)*

**Burrows–Wheeler transform (BWT) / FM-index** — a compression transform
repurposed as a search index; makes exact-match lookup O(query length),
independent of genome size. The engine inside `bwa`. *(Ch. 7)*

**Central dogma** — DNA → RNA → protein; sequence information flows
forward. *(Ch. 2)*

**CIGAR** — run-length-encoded edit script describing how a read aligns
to the reference (`50M2D48M`); `M` means aligned, not matching. *(Ch. 9)*

**Codon** — three consecutive RNA bases, mapping via the genetic code to
one amino acid or stop. *(Ch. 2)*

**Consensus / profile matrix** — per-column base counts across aligned
sequences; the counting form of a position weight matrix. *(Ch. 5)*

**Coverage / depth** — how many reads overlap a position, on average
("30×"); the denominator of confidence in every call. *(Ch. 7, 8, 14)*

**Degeneracy** — 64 codons → 21 outcomes: multiple codons per amino acid,
usually differing in the third position; why synonymous variants exist.
*(Ch. 2)*

**Exon / intron** — in eukaryotic genes, the retained coding chunks
(exons) and the intervening sequence (introns) removed by splicing.
*(Ch. 3)*

**F1** — harmonic mean of precision and recall. *(Ch. 15)*

**FASTA** — sequence format: `>` header plus wrapped sequence lines.
*(Ch. 5, 8)*

**FASTQ** — read format: four lines per read — id, bases, `+`, per-base
Phred quality string. *(Ch. 8)*

**FLAG** — the SAM bitfield of per-read booleans (reverse strand,
unmapped, duplicate, pairing…). *(Ch. 9)*

**Frameshift** — an indel of length not divisible by 3 inside coding
sequence; shifts every downstream codon boundary. *(Ch. 4)*

**Gene** — a stretch of DNA that gets transcribed; in humans, interrupted
by introns and often multiply spliced. *(Ch. 3)*

**Genome build** — a versioned release of the reference genome
(hg19/GRCh37, hg38/GRCh38); coordinates are meaningless without one.
*(Ch. 11)*

**Genotype (GT)** — which alleles a sample carries at a site: `0/0`,
`0/1`, `1/1`; `|` when phased. *(Ch. 10)*

**GFF / GTF** — rich annotation formats (gene models: exons →
transcripts → genes, with strand and frame); **1-based closed**
coordinates. *(Ch. 11, 12)*

**GIAB (Genome in a Bottle)** — NIST-led consortium producing benchmark
genomes (HG002 et al.): truth VCFs plus high-confidence BEDs. *(Ch. 13)*

**GWAS** — genome-wide association study: regress a trait on each variant
across the genome; genome-wide significance p < 5×10⁻⁸; PCs included as
covariates against stratification. *(Ch. 17)*

**Hamming distance (HAMM)** — count of mismatching positions between
equal-length strings; shattered by a single indel. *(Ch. 5)*

**Haplotype** — the sequence actually carried on one chromosome copy;
haplotype comparison is the correct way to compare VCFs. *(Ch. 8, 15)*

**Heterozygous / homozygous** — the two chromosome copies differ / agree
at a site. *(Ch. 4)*

**HGVS** — the standard variant nomenclature (`g.`, `c.`, `p.`); protein
residues numbered from the initiator Met = 1 (vs mature numbering, one
lower). *(Ch. 4)*

**High-confidence regions** — the ~90% of the genome where GIAB claims
truth; calls outside are unscoreable, not wrong. *(Ch. 13)*

**Indel** — small insertion or deletion (up to ~50bp). *(Ch. 4)*

**LD (linkage disequilibrium)** — correlation between nearby variants,
inherited in blocks; pruned before PCA (it's multicollinearity).
*(Ch. 16)*

**MAF (minor allele frequency)** — population frequency of a site's rarer
allele; MAF > 5% = "common variant". *(Ch. 16)*

**MANE Select** — agreed canonical transcript per gene (NCBI + EMBL-EBI).
*(Ch. 4)*

**MAPQ** — Phred-scaled confidence that a read's *placement* is correct;
0 = maps equally well elsewhere (repeats). *(Ch. 9)*

**Motif** — a short recurring pattern (e.g. a binding site); modelled
probabilistically via PWMs. *(Ch. 5)*

**Needleman–Wunsch / Smith–Waterman** — global / local alignment by
dynamic programming. *(Ch. 6)*

**ORF (open reading frame)** — in some frame: AUG, a run of codons with
no stop, then a stop; long ORFs are statistical evidence of a gene.
*(Ch. 3)*

**PCA (in popgen)** — eigendecomposition of the standardised genotype
matrix; top components recover ancestry/geography. *(Ch. 16, 17)*

**Phred score** — Q = −10·log₁₀ P(error), stored as ASCII with offset 33;
used for base quality, MAPQ, and VCF QUAL alike. *(Ch. 8)*

**Pileup** — all read bases stacked over one reference position; the
variant caller's input. *(Ch. 10)*

**Precision / recall** — TP/(TP+FP) and TP/(TP+FN): "how much of what I
called was real" / "how much of what was real did I find". *(Ch. 15)*

**Read** — one sequenced fragment (~150bp for standard Illumina); paired
ends = both ends of one ~350bp fragment. *(Ch. 8, 14)*

**Reading frame** — one of six ways (3 offsets × 2 strands) to partition
a sequence into codons. *(Ch. 3)*

**Reference genome** — the agreed consensus sequence everything is
reported against; a mosaic of donors, not a person. *(Ch. 4)*

**Reverse complement** — the partner strand written 5′→3′: complement the
bases, reverse the string; the most-used operation in the field. *(Ch. 1)*

**rs334** — the sickle-cell variant: hg38 chr11:5,227,002 T>A (genome
plus strand), c.20A>T / p.Glu7Val on HBB; this book's running example.
*(Ch. 4, 8, 10, 15)*

**SAM** — text format for alignments: 11 mandatory columns per read.
*(Ch. 9)*

**SNV / SNP** — single-base variant; SNP historically implies "common in
a population". *(Ch. 4)*

**Soft clip (S in CIGAR)** — read bases present but unaligned at a read's
edge; piles of them mark structural-variant breakpoints. *(Ch. 9)*

**Splicing** — cutting introns out of the transcript and joining exons
before translation. *(Ch. 3)*

**Strand (+ / −)** — the two antiparallel DNA strands; features and reads
belong to one; sequences are always written 5′→3′. *(Ch. 1)*

**Structural variant (SV)** — large (>50bp) rearrangement: deletion,
duplication, inversion, translocation. *(Ch. 4)*

**Synonymous / missense / nonsense / start-lost / stop-lost** — the
consequence taxonomy of coding SNVs. *(Ch. 4)*

**Transcription / translation** — DNA → RNA (string-level: T→U) / RNA →
protein (codon table, three at a time). *(Ch. 2)*

**Trio** — child plus both parents; inheritance consistency as a truth
signal (HG002/3/4). *(Ch. 13)*

**Variant representation problem** — the same physical event (especially
indels in repeats) legally written as different VCF records; addressed by
normalisation and haplotype comparison. *(Ch. 15)*

**VCF** — variant format: one record per site where a sample differs from
the reference; a *diff*, meaningless without its build; absence means
"reference", not "no data". *(Ch. 10)*
