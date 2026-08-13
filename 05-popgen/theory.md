# 05 — Population structure falls out of a PCA

The module where your existing skills stop being *transferable* and simply
*are* the skill: everything is linear algebra on a matrix, and the biology
is in the interpretation.

The theory now lives in **Part V of the book**
([phreakv6.github.io/genomics](https://phreakv6.github.io/genomics/)):

- [Ch. 16 — Two Thousand Genomes as a Matrix](../book/src/ch16-genomes-as-a-matrix.md) —
  1000 Genomes, the 0/1/2 encoding, MAF filtering, LD pruning as
  multicollinearity, Patterson normalisation
- [Ch. 17 — Geography in the Genome](../book/src/ch17-geography-in-the-genome.md) —
  why PCA recovers geography, reading the plot (including the within-SAS
  cline), the caveats, and the GWAS/factor-model connection

## Then

`fetch_data.sh` slices real 1000 Genomes chromosome 22 data out of the
remote VCF; `pca.py` does the filtering, pruning, and standardisation with
`scikit-allel` and plots the result coloured by population.

```bash
bash 05-popgen/fetch_data.sh
python 05-popgen/pca.py
```
