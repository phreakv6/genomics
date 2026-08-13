# Chapter 12 — Annotation: GFF and GTF *(parked)*

> **Status: deliberately unwritten.** This chapter is a placeholder — the
> outline of work the repo has planned but not yet done. It's included
> rather than hidden because the gap it marks is real and load-bearing,
> and pretending the book covers annotation would violate its own rules.

## Why this chapter needs to exist

Annotation files are where three of this book's recurring characters
actually live:

- **Strand.** Chapter 8 noted that HBB is transcribed from the minus
  strand — which is why the genome says `T>A` where the coding sequence
  says `A>T`. That fact was *asserted in a code comment*. The place it's
  actually recorded, authoritatively, is the strand column of an
  annotation file.
- **Exon structure.** Chapter 3 explained that human genes are interrupted
  and that naive ORF finding therefore fails. The exon/intron boundaries —
  the information that repairs it — are annotation.
- **Transcripts and frames.** Chapter 4's deepest point was that a
  variant's consequence depends on the annotation (which transcript, which
  frame). GFF/GTF is the file format that plurality of knowledge is
  shipped in — one line per feature, with parent/child relationships
  linking exons → transcripts → genes.

The repo currently touches none of this in code: modules 01–05 use
sequence, alignments, and variants, but never read an annotation file.
Every strand and frame fact in the book so far was hardcoded from
knowledge. That's the gap.

## The planned demonstration

When this chapter gets written, its spine is already designed — a
round-trip proof in the spirit of Chapter 3's ORF check:

1. Pull the real HBB annotation (GTF) from Ensembl.
2. Read the **strand column** — turning the "HBB is minus-strand"
  assertion into data.
3. Rebuild the coding sequence from the **exon list**: fetch each exon's
  genomic sequence, stitch in transcript order, reverse-complement
  (minus strand), and trim to the CDS.
4. Show the result **matches the `HBB_CDS` constant** hardcoded in module
  01 since Chapter 1 — closing the loop between the book's first hardcoded
  string and the file format that makes it derivable.

## Until then

The three-sentence version of the format, from Chapter 10: GFF/GTF is the
rich interval format — tab-separated features with chromosome, source,
type (gene/transcript/exon/CDS), **1-based closed** coordinates
(Chapter 11!), score, **strand**, frame, and a semi-structured attributes
column linking features into gene models. `bedtools` speaks it;
VEP/SnpEff consume it to annotate variants; Ensembl and GENCODE publish
the canonical human ones.

## Further reading

- The [Ensembl GFF3/GTF file documentation](https://www.ensembl.org/info/website/upload/gff.html)
  — short, and enough to read a gene model by hand.
