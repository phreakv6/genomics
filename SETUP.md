# Setup

Everything here runs on a laptop. No cluster, no wet lab, no institutional
access, no paid data. Total download across all five modules is ~55MB; total
disk once indexed is ~500MB.

## 1. Get conda

If you don't have it, [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
or [Miniforge](https://github.com/conda-forge/miniforge) — either is fine.
Miniforge if you'd rather avoid Anaconda's licensing entirely.

Conda rather than pip because most of these are compiled C programs, not Python
packages. `bioconda` is how the entire field distributes them, and fighting that
convention is not a hill worth choosing.

## 2. Create the environment

**Apple Silicon Mac (M1/M2/M3/M4) — read this bit:**

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate bio
conda config --env --set subdir osx-64     # so later installs stay consistent
```

Much of bioconda has no `osx-arm64` build. Without `CONDA_SUBDIR=osx-64`, conda
reports the packages as simply not existing — a confusing error for a package
that obviously exists. The Intel builds run under Rosetta. The first Python
import takes ~50 seconds while translations are cached; after that it's normal
speed.

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

Every `run.sh` checks its tools are present and tells you what to do if they
aren't, so you can't get far with a broken environment without noticing.

## 3. Run the modules

Modules 01 and 02 need only Python — no environment, no downloads. Start there.

```bash
# 01 -- central dogma, from scratch
python 01-central-dogma/central_dogma.py

# 02 -- Rosalind problems, then alignment by dynamic programming
python 02-rosalind/rosalind.py
python 02-rosalind/alignment.py

# 03 -- formats.  ~4KB download, a few seconds
python 03-formats/make_data.py     # fetches the HBB locus, simulates reads
bash   03-formats/run.sh           # FASTQ -> BAM -> VCF
python 03-formats/dissect.py       # take every format apart
python 03-formats/coordinates.py   # the off-by-one footgun, concretely

# 04 -- real GIAB data.  ~60MB download, ~460MB on disk, ~4 min total
bash   04-pipeline/fetch_data.sh   # chr11 reference + real HG002 reads
bash   04-pipeline/run.sh          # align, call, compare to truth
python 04-pipeline/benchmark.py    # precision/recall + error analysis

# 05 -- population structure.  ~10MB download, ~2 min
bash   05-popgen/fetch_data.sh     # 1000 Genomes slice
python 05-popgen/pca.py            # PCA, writes data/pca.png
```

Read each module's `theory.md` before running its code. The code is written to
be read alongside the prose — the functions map onto the sections.

## Costs, so nothing surprises you

| Module | Download | Disk | Time |
|--------|----------|------|------|
| 01, 02 | none | none | instant |
| 03 | 4KB | 1MB | ~10s |
| 04 | ~60MB | ~460MB | ~4 min (90s of it is `bwa index`) |
| 05 | ~10MB | ~12MB | ~2 min |

All downloaded data lands in each module's `data/` directory, which is
gitignored. Delete those directories any time; the fetch scripts are re-runnable
and skip work that's already done.

## Data sources, all public

- **UCSC Genome Browser** — reference sequence (hg38) and its REST API
- **Genome in a Bottle (NIST)** — HG002 reads and the benchmark truth set
- **1000 Genomes Project** — Phase 3 genotypes for 2,504 people

Nothing here requires an account. Modules 04 and 05 stream slices out of huge
remote files using the index, so you never download the 250GB BAM or the 200MB
chromosome VCF.

## Troubleshooting

**`PackagesNotFoundError` on an Apple Silicon Mac** — you forgot
`CONDA_SUBDIR=osx-64`. Delete the environment and recreate it.

**A download hangs or fails** — NCBI and EBI's public FTP endpoints throttle
under load. The fetch scripts are re-runnable and skip completed steps; just run
again.

**Zero results, no error** — the classic genomics failure. Suspect chromosome
naming (`chr11` vs `11`) or a reference build mismatch before you suspect
biology. See footgun #2 in `03-formats/theory.md`.

**`bwa index` seems stuck** — it's ~90 seconds of silence for chr11. Normal.

## If you want to learn along

Issues and PRs welcome. The most useful contributions are corrections — this is
a repo written while learning, so if something is wrong or misleading, say so.
