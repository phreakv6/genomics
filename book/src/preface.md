# Preface — How to Read This Book

This book is the theory behind [`phreakv6/genomics`](https://github.com/phreakv6/genomics),
a repo built while learning genomics from scratch: five modules that run from
DNA-as-a-string all the way to a real variant-calling pipeline scored against a
truth set, and a PCA that recovers human geography from raw genotypes. The repo
holds the runnable code; this book holds everything you need to know to
understand every line of it.

## Who it's for

Someone who can already program and knows statistics, and knows nothing about
biology. That is exactly the position the author was in, and the book is
written to that reader:

- **Nothing about programming is explained.** If you need Python explained,
  this is the wrong book. Dictionaries, dynamic programming, and bit flags
  appear without ceremony.
- **Everything about biology is explained**, from "what is DNA" upward, the
  first time it appears. No prior vocabulary is assumed — when a term like
  *codon*, *ORF*, or *Phred score* first shows up, it gets defined properly,
  not gestured at.
- **The *why* matters as much as the *what*.** The field does many things that
  look arbitrary — 1-based coordinates here, 0-based there, quality scores as
  ASCII characters, alignment via a compression algorithm's transform. Where a
  convention exists, this book tries to say why, because conventions you
  understand are conventions you stop tripping over.

## How the book and the code fit together

Every chapter is backed by code in the repo, and the two are kept deliberately
in sync:

- **Code blocks are quoted from the repo's actual source** — the same
  functions you can run, not simplified pedagogical forks.
- **Output blocks are captured from real runs.** When the book says a 93bp
  sequence gives frames of 31, 30, and 30 residues, that number came out of
  `python central_dogma.py`, not out of anyone's head. In a teaching text the
  reader trusts the prose over the code, so prose that drifts from what the
  code does is the most expensive possible failure — this book's rule is that
  numbers are computed, never asserted.
- **Each chapter ends with a "Run it" section**: the exact commands that
  produce everything the chapter showed you.

The five parts map onto the five modules:

| Part | Chapters | Module | What you'll be able to do afterwards |
|------|----------|--------|--------------------------------------|
| I — The Central Dogma | 1–4 | `01-central-dogma/` | Read DNA → RNA → protein as string operations; classify a mutation |
| II — Comparing Sequences | 5–7 | `02-rosalind/` | Implement and reason about alignment; know why `bwa` exists |
| III — The Formats | 8–12 | `03-formats/` | Read FASTA/FASTQ/SAM/BAM/VCF fluently; dodge the coordinate traps |
| IV — A Real Pipeline | 13–15 | `04-pipeline/` | Call variants on real human data and score yourself honestly |
| V — Populations | 16–17 | `05-popgen/` | Turn 2,504 genomes into a matrix and interpret its principal components |

Parts I and II need nothing but Python. Parts III–V use a conda environment
and public data — about 55MB of downloads total, all on a laptop. The
[Setup chapter](setup.md) covers it.

## Ground rules, inherited from the repo

- **Public data only.** 1000 Genomes, Genome in a Bottle, UCSC. Nothing needs
  an account; nothing is anyone's private genome.
- **Every result is checked against a truth set where one exists.** The reason
  Part IV uses the GIAB benchmark sample is precisely that you can score
  yourself against it — and the book reports the actual scores, including the
  errors, because the errors turn out to be the most instructive part.
- **Corrections are the most valuable contribution.** This is a book written
  while learning. If something is wrong, [open an issue](https://github.com/phreakv6/genomics/issues).

## A note on reading order

The chapters build strictly forward — later parts lean on earlier vocabulary
without re-defining it. If you already know some biology, Part I will be fast
but is still worth skimming for the conventions (HGVS numbering, coordinate
systems, strand notation) that the rest of the book uses without comment.
There is a [Glossary](glossary.md) at the back for when a term has scrolled
out of memory.
