# Setup — The Environment

Everything in this book runs on a laptop. No cluster, no wet lab, no
institutional access, no paid data. Total download across all five parts is
~55MB; total disk once indexed is ~500MB.

Parts I and II (chapters 1–7) need **nothing but Python** — you can read the
next seven chapters and run every example with a stock interpreter. Come back
here when you reach Part III.

## Get conda

If you don't have it: [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
or [Miniforge](https://github.com/conda-forge/miniforge) — either is fine;
Miniforge if you'd rather avoid Anaconda's licensing entirely.

Why conda rather than pip: most genomics tools (`samtools`, `bcftools`, `bwa`)
are compiled C programs, not Python packages. [Bioconda](https://bioconda.github.io/)
is how the entire field distributes them, and fighting that convention is not
a hill worth choosing.

## Create the environment

**Apple Silicon Mac (M1/M2/M3/M4) — read this bit:**

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate bio
conda config --env --set subdir osx-64     # so later installs stay consistent
```

Much of bioconda has no `osx-arm64` build. Without `CONDA_SUBDIR=osx-64`,
conda reports the packages as simply not existing — a confusing error for a
package that obviously exists. The Intel builds run under Rosetta. The first
Python import takes ~50 seconds while translations are cached; after that it's
normal speed.

**Linux, Intel Mac, or WSL:**

```bash
conda env create -f environment.yml
conda activate bio
```

**Verify:**

```bash
samtools --version | head -1
python -c "import allel, pysam; print(allel.__version__, pysam.__version__)"
```

Every `run.sh` in the repo checks its tools are present and tells you what to
do if they aren't, so you can't get far with a broken environment without
noticing.

## Costs, so nothing surprises you

| Part | Module | Download | Disk | Time |
|------|--------|----------|------|------|
| I, II | 01, 02 | none | none | instant |
| III | 03 | 4KB | 1MB | ~10s |
| IV | 04 | ~60MB | ~460MB | ~4 min (90s of it is `bwa index`) |
| V | 05 | ~10MB | ~12MB | ~2 min |

All downloaded data lands in each module's `data/` directory, which is
gitignored. Delete those directories any time; the fetch scripts are
re-runnable and skip work that's already done.

## Data sources, all public

- **UCSC Genome Browser** — reference sequence (hg38) and its REST API
- **Genome in a Bottle (NIST)** — HG002 reads and the benchmark truth set
- **1000 Genomes Project** — Phase 3 genotypes for 2,504 people

Nothing requires an account. Parts IV and V stream slices out of huge remote
files using their indexes, so you never download the 250GB BAM or the 200MB
chromosome VCF — how that trick works is explained in
[Chapter 14](ch14-the-pipeline.md).

## Troubleshooting

**`PackagesNotFoundError` on an Apple Silicon Mac** — you forgot
`CONDA_SUBDIR=osx-64`. Delete the environment and recreate it.

**A download hangs or fails** — NCBI and EBI's public FTP endpoints throttle
under load. The fetch scripts are re-runnable and skip completed steps; just
run again.

**Zero results, no error** — the classic genomics failure. Suspect chromosome
naming (`chr11` vs `11`) or a reference build mismatch before you suspect
biology. [Chapter 11](ch11-coordinates-and-builds.md) is about exactly this.

**`bwa index` seems stuck** — it's ~90 seconds of silence for chr11. Normal.
