# Setup

Everything here runs on a laptop: no cluster, no wet lab, no accounts, no
paid data. ~55MB of downloads, ~500MB of disk, about 10 minutes.

The full setup guide — conda, the **Apple Silicon gotcha**, per-module
costs, data sources, troubleshooting — is the book's Setup chapter:

**→ [Setup — The Environment](https://phreakv6.github.io/genomics/setup.html)**
(source: [`book/src/setup.md`](book/src/setup.md), renders on GitHub)

The short version:

```bash
# Apple Silicon Mac (M1-M4) -- the CONDA_SUBDIR prefix is required:
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate bio
conda config --env --set subdir osx-64

# Linux / Intel Mac / WSL:
conda env create -f environment.yml
conda activate bio
```

Modules 01–02 need no environment at all — just Python. Start there.

If something breaks, the Setup chapter's troubleshooting section covers the
four failures people actually hit (`PackagesNotFoundError`, hanging
downloads, silent zero-result runs, `bwa index` "stuck").
