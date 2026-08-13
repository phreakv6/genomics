# Rosalind Problem Index

The classic [Rosalind](https://rosalind.info) problems this book solves,
where each is implemented in the repo, and which chapter explains it.

| Code | Problem | Function | Chapter |
|------|---------|----------|---------|
| `DNA` | Count nucleotides | `central_dogma.py::count_bases` | [Ch. 1](ch01-dna-is-a-string.md) |
| `RNA` | Transcribe DNA → RNA | `central_dogma.py::transcribe` | [Ch. 2](ch02-dna-to-protein.md) |
| `REVC` | Reverse complement | `central_dogma.py::reverse_complement` | [Ch. 1](ch01-dna-is-a-string.md) |
| `GC` | Highest GC content in a FASTA | `central_dogma.py::gc_content` + `rosalind.py::parse_fasta` | [Ch. 5](ch05-distance-before-alignment.md) |
| `HAMM` | Hamming distance | `rosalind.py::hamming` | [Ch. 5](ch05-distance-before-alignment.md) |
| `SUBS` | Motif occurrences (overlapping) | `rosalind.py::find_motif` | [Ch. 5](ch05-distance-before-alignment.md) |
| `PROT` | Translate RNA → protein | `central_dogma.py::translate` | [Ch. 2](ch02-dna-to-protein.md) |
| `CONS` | Consensus & profile matrix | `rosalind.py::profile_matrix`, `consensus` | [Ch. 5](ch05-distance-before-alignment.md) |
| `ORF`* | Open reading frames | `central_dogma.py::find_orfs` | [Ch. 3](ch03-reading-frames-and-orfs.md) |

\* The repo solves the substance of `ORF` (six-frame ORF finding on both
strands) without having been written against Rosalind's exact input/output
format.

Alignment (Chapter 6) corresponds to Rosalind's later `EDIT`/`GLOB`/`LOCA`
problems; the repo implements Needleman–Wunsch and Smith–Waterman directly
in `alignment.py` rather than through Rosalind's framing.

If you want more practice, the natural next Rosalind problems after this
book's set: `PERM`-tier warmups aside, try `LCSM` (shared motifs), `TRAN`
(transition/transversion ratio — a nice QC statistic), and `LONG`
(fragment assembly, a taste of genome assembly proper).
