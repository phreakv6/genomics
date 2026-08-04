# Working on this repo

## What this is

A learning repo, not a product. Bharani is a software engineer and investor of
10+ years learning genomics from scratch, module by module, with Claude Code.
Public data only, laptop only. It is on GitHub (`phreakv6/genomics`, public) so
others can learn alongside.

The five modules run central dogma → Rosalind/alignment → file formats → a real
GIAB variant-calling pipeline → 1000 Genomes PCA. Each has a `theory.md` and
runnable code. See `README.md` for results and `SETUP.md` for environment.

## How he wants to work

**Currently working through module 01.** Do not rush ahead. He asked
explicitly not to overcomplicate things while still early in the sequence.

- **Assume strong engineering, zero biology.** He does not need Python, stats,
  data structures, or shell explained. He does need biology vocabulary and
  domain convention explained, and he appreciates knowing *why* the field does
  something the way it does, not just what it does.
- **Answer conceptual questions with a runnable demo.** When he asks "what
  is X" or "wouldn't Y break", the good answer computes the thing and shows
  real output, then explains. Several of the best moments in this repo came
  from exactly that.
- **He is here to understand, not to type.** He said up front: the point is not
  to write the code manually, it is to understand every line and why it exists.
  So write the code, run it, and walk through the output.
- Analogies to investing/markets land well (LD pruning as multicollinearity,
  GIAB high-confidence regions as survivorship bias, PCs as factor models).

## The rule that matters most here

**Never let prose drift from what the code does.** His questions have caught
three real errors so far, and every one was a comment or a display lying about
correct code:

- frame output clipped to `[:30]`, hiding that frames have unequal lengths
- `p.Glu6Val` labelled HGVS when HGVS is `p.Glu7Val` (mature vs HGVS numbering)
- `# 12 aa + stop` when it is 12 codons, i.e. 11 aa + stop

In a teaching repo a learner trusts the prose over the code, so this is the
expensive failure mode. Prefer computing numbers into the output over asserting
them in comments where they rot silently. Verify claims by running them.

## Environment

```
conda activate bio          # built with CONDA_SUBDIR=osx-64 (Apple Silicon)
```

Much of bioconda has no `osx-arm64` build. First Python import takes ~50s while
Rosetta caches; after that it is fast. Modules 01–02 need only plain Python.
`data/` dirs are gitignored downloads, safe to delete, re-fetchable.

Note: `bcftools`/`samtools` cache remote `.tbi`/`.bai` index files into the
working directory when doing remote indexed access. They are gitignored, but
that is why stray index files appear in the repo root.

## Conventions

- Reference build goes in every filename: `sample.hg38.bam`.
- Coordinates: state the convention whenever it could be ambiguous. Module 03
  documents the 0-based/1-based zoo and the hg19/hg38 trap.
- Commit messages: explain what was wrong and why it mattered, not just what
  changed. Push to `main` when he asks; he does the asking.

## Parked / open threads

- **Annotation (GFF/GTF) is missing entirely from the repo.** Modules 03–05 use
  sequence, alignments and variants but never an annotation file — which is
  where strand, exons and gene names actually live. Planned addition to module
  03: pull the real HBB GTF, read the strand column, rebuild the CDS from the
  exon list, and show it matching the `HBB_CDS` constant hardcoded in module 01.
  **Deliberately parked** until he is further along.
- `03-formats/dissect.py` closing note still says `p.Glu6Val`, inconsistent
  with module 01's corrected HGVS numbering. Small, unfixed.
- `check_env.sh`'s Linux branch and the non-`CONDA_SUBDIR` setup path are
  written from knowledge, never tested. First Linux learner is the real test.

## Origin

The project came from a conversation about where a multi-disciplinary
generalist could find edge in an underexplored field. The conclusion was that
edge is a *shape* — under-covered, messy data, cross-domain synthesis,
build-your-own-tools, no need to publish — not a glamorous field. Genomics
tooling was picked as the buildable on-ramp because the data is public and the
gap between "published model" and "thing a scientist can use" is real. Longer
term the higher-edge play is the pharma-intelligence × biology seam he already
occupies professionally.
