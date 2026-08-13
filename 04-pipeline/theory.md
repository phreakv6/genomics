# 04 — The real pipeline, and scoring yourself against truth

Module 03 was a rehearsal. Here the same pipeline runs on **real reads
from a real person** (GIAB's HG002), and then — the point of the module —
we **measure how wrong we were**.

The theory now lives in **Part IV of the book**
([phreakv6.github.io/genomics](https://phreakv6.github.io/genomics/)):

- [Ch. 13 — What Truth Looks Like: HG002 and GIAB](../book/src/ch13-hg002-and-giab.md) —
  how a truth set is made; why the high-confidence BED is the
  intellectually important file (survivorship bias)
- [Ch. 14 — FASTQ → BAM → VCF, For Real](../book/src/ch14-the-pipeline.md) —
  every pipeline step with its reason; the remote-slicing trick
- [Ch. 15 — Keeping Score](../book/src/ch15-keeping-score.md) —
  precision/recall/F1, the results, every error read individually, the
  variant-representation problem

## Then

`fetch_data.sh` pulls chr11 and slices real HG002 reads out of a remote
250GB BAM; `run.sh` realigns and calls them; `benchmark.py` scores the
result against the GIAB truth set restricted to high-confidence regions
and inspects every error.

```bash
bash 04-pipeline/fetch_data.sh && bash 04-pipeline/run.sh
python 04-pipeline/benchmark.py
```
