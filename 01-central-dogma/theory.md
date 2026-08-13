# 01 — The Central Dogma

The whole of module 1 is one sentence: **DNA → RNA → protein**.

The theory for this module now lives in **Part I of the book**
([phreakv6.github.io/genomics](https://phreakv6.github.io/genomics/)),
which covers it at full depth with figures and the real output of this
module's code:

- [Ch. 1 — DNA Is a String (Almost)](../book/src/ch01-dna-is-a-string.md) —
  bases, pairing, 5′→3′ polarity, reverse complement, GC content
- [Ch. 2 — From DNA to Protein](../book/src/ch02-dna-to-protein.md) —
  transcription, the genetic code, translation
- [Ch. 3 — Reading Frames and ORFs](../book/src/ch03-reading-frames-and-orfs.md) —
  six frames, ORF finding, why it fails on real genomes
- [Ch. 4 — When One Base Changes](../book/src/ch04-when-one-base-changes.md) —
  variant taxonomy, sickle cell, HGVS numbering, frames as annotation

## Then

`central_dogma.py` implements every operation from those chapters from
scratch — complement, reverse complement, transcription, the codon table,
translation, six-frame translation, ORF finding — and runs them on a real
human gene fragment. Each function maps to a book section.

```bash
python 01-central-dogma/central_dogma.py
```
