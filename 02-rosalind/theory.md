# 02 — Rosalind, and the algorithm underneath everything

[rosalind.info](https://rosalind.info) is Project Euler for bioinformatics. The first dozen problems
walk you through the central dogma as code — you already have that from module 01. What's worth
your time is the *second* thing in this module: sequence alignment, which is the algorithmic core of
the entire field and the thing module 04 will be running at scale.

## The classic problems (and what each is really about)

| Code | Problem | The real point |
|------|---------|----------------|
| `DNA` | Count nucleotides | Hello world |
| `RNA` | Transcribe | T→U |
| `REVC` | Reverse complement | The most-used operation in the field |
| `GC` | Highest GC content in a FASTA | Your first parser. FASTA is everywhere |
| `HAMM` | Hamming distance | Counting mismatches = the crudest possible alignment score |
| `SUBS` | Find motif occurrences | Overlapping matches — `str.find` in a loop, not regex |
| `PROT` | Translate RNA | Codon table |
| `CONS` | Consensus & profile matrix | Position weight matrices — the basis of motif finding |

`HAMM` is more load-bearing than it looks. It only works on equal-length strings with no gaps. The
moment you allow insertions and deletions — which biology does constantly — you need alignment, and
alignment needs dynamic programming.

## Alignment: the actual core

**The problem.** Given two sequences, find the correspondence between them that maximises a score,
where you get points for matches, lose points for mismatches, and lose points for gaps (insertions
or deletions). Biologically: two sequences descend from a common ancestor, and you're reconstructing
the cheapest set of mutations that turns one into the other.

**Needleman–Wunsch (1970) — global alignment.** Align the sequences end to end. Build an
(m+1)×(n+1) matrix where cell `(i,j)` holds the best score for aligning the first `i` characters of
one against the first `j` of the other. Each cell is the best of three choices:

```
F(i,j) = max( F(i-1, j-1) + score(a[i], b[j])   # diagonal: align the two characters
              F(i-1, j)   + gap                 # up:   consume from a, gap in b
              F(i,   j-1) + gap )               # left: consume from b, gap in a
```

Then trace back from the bottom-right corner to recover the alignment itself. This is the same
edit-distance DP you already know from computer science, with a biology-flavoured scoring function.
O(mn) time and space.

**Smith–Waterman (1981) — local alignment.** One change: clamp negative scores to zero, and start
the traceback from the highest cell anywhere in the matrix rather than the corner. That makes it
find the best-matching *sub*region instead of forcing an end-to-end alignment. This is what you
want when a short read overlaps part of a long genome, or when two proteins share one conserved
domain amid otherwise unrelated sequence.

**Scoring matters.** For DNA, a simple +1/−1 works. For proteins it doesn't: substituting leucine
for isoleucine (both greasy, similar size) is far more forgivable than leucine for aspartate.
So protein alignment uses substitution matrices — BLOSUM62, PAM — whose entries are log-odds derived
from observed substitution frequencies in real alignments. Positive = happens more often than
chance. It's an empirical prior, not a theory.

Also: real gaps are not linear. Opening a gap is rare; extending an existing one is cheap, because
one mutational event deletes many bases at once. Hence **affine gap penalties** (`open` + `k ×
extend`), which is what real tools use.

## Why you can't use this on a genome

Aligning one 100bp read against the 3.1-billion-base human genome with Smith–Waterman costs ~3×10¹¹
cell computations. Now do that for 500 million reads. It's off by many orders of magnitude.

The fix, and the reason modern aligners exist:

- **Seed-and-extend** (BLAST, 1990) — find short exact matches first, then do expensive DP only
  around those seeds. Heuristic; sacrifices the optimality guarantee for ~50× speed.
- **Burrows–Wheeler transform + FM-index** (BWA, Bowtie, ~2009) — this is the real breakthrough.
  Build a compressed index of the genome once, and exact-match queries become O(query length),
  independent of genome size. The whole human genome index fits in a few GB of RAM. This is why
  `bwa` in module 04 aligns millions of reads on your laptop. The BWT is the same transform used in
  `bzip2` — a compression algorithm turned into a search index. Worth appreciating: the single
  biggest speedup in genomics came from a data-structures trick, not from biology.

So the mental model going into module 04: `bwa` uses BWT/FM-index to *find candidate locations*
fast, then runs Smith–Waterman-style DP only in those small windows to get the precise alignment
with gaps. The code in this module is the second half of that sentence, written out.

## What to read

- Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach* — the DP chapters.
  This book teaches by making you implement; it pairs with Rosalind directly.
- Durbin, Eddy, Krogh & Mitchison, *Biological Sequence Analysis* — Ch. 2 (alignment) and the HMM
  chapters. Mathematically serious; a second-pass book, but you'll enjoy it.
- Buffalo, *Bioinformatics Data Skills* — start it now, alongside module 03.

## Then

`rosalind.py` solves the classic problems; `alignment.py` implements Needleman–Wunsch and
Smith–Waterman with traceback and prints the matrices small enough to actually read.
