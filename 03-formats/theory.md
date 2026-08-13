# 03 — Formats, coordinates, and the footguns

The module with no intellectual glamour and 90% of real bugs.

The theory now lives in **Part III of the book**
([phreakv6.github.io/genomics](https://phreakv6.github.io/genomics/)):

- [Ch. 8 — Sequences on Disk: FASTA and FASTQ](../book/src/ch08-fasta-and-fastq.md) —
  the formats, Phred scores, why variant calling is Bayesian
- [Ch. 9 — Alignments on Disk: SAM and BAM](../book/src/ch09-sam-and-bam.md) —
  the 11 columns, FLAG, MAPQ, CIGAR, sort-and-index
- [Ch. 10 — Variants on Disk: VCF](../book/src/ch10-vcf.md) —
  the pileup, record anatomy, "a diff against the reference"
- [Ch. 11 — Coordinates, Builds, and Silent Wrongness](../book/src/ch11-coordinates-and-builds.md) —
  0- vs 1-based, hg19 vs hg38, `chr11` vs `11`, the five habits
- [Ch. 12 — Annotation: GFF/GTF](../book/src/ch12-annotation.md) —
  parked; the planned HBB-GTF demonstration lives there as an outline

## Then

`make_data.py` fetches the real hg38 HBB locus and simulates reads with
rs334 planted; `run.sh` runs the canonical FASTQ → BAM → VCF six commands;
`dissect.py` takes every format apart field by field; `coordinates.py`
makes the off-by-one footgun concrete.

```bash
python 03-formats/make_data.py && bash 03-formats/run.sh
python 03-formats/dissect.py && python 03-formats/coordinates.py
```
