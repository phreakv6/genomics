# 02 — Rosalind, and the algorithm underneath everything

[rosalind.info](https://rosalind.info) is Project Euler for
bioinformatics; the real content of this module is **sequence alignment**,
the algorithmic core of the field.

The theory now lives in **Part II of the book**
([phreakv6.github.io/genomics](https://phreakv6.github.io/genomics/)):

- [Ch. 5 — Distance Before Alignment](../book/src/ch05-distance-before-alignment.md) —
  FASTA parsing (GC), Hamming distance (HAMM), overlapping motifs (SUBS),
  profile matrices (CONS), and why one indel breaks all of it
- [Ch. 6 — Alignment](../book/src/ch06-alignment.md) —
  Needleman–Wunsch and Smith–Waterman, scoring, BLOSUM, affine gaps
- [Ch. 7 — Why You Can't Align a Genome This Way](../book/src/ch07-scaling-alignment.md) —
  the scaling arithmetic, seed-and-extend, BWT/FM-index; the `bwa` mental model

## Then

`rosalind.py` solves the classic problems; `alignment.py` implements
Needleman–Wunsch and Smith–Waterman with traceback and prints the matrices
small enough to actually read.

```bash
python 02-rosalind/rosalind.py
python 02-rosalind/alignment.py
```
