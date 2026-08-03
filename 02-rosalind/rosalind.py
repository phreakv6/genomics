"""The classic Rosalind problems, with real data rather than toy strings.

Run: python 02-rosalind/rosalind.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "01-central-dogma"))
from central_dogma import transcribe, translate, reverse_complement, gc_content  # noqa: E402


# ---------------------------------------------------------------------------
# FASTA -- the format you will parse a thousand times
# ---------------------------------------------------------------------------
# >header line describing the sequence
# ACGTACGTACGT...      <- sequence, wrapped at 60 or 70 chars
# ACGTACGT
# >next_header
# ...
#
# That's it. The wrapping is why you cannot just read line-by-line and treat
# each line as a record -- you must accumulate until the next '>'.

def parse_fasta(text):
    """Yield (header, sequence) pairs. Handles multi-line sequences."""
    header, chunks = None, []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


SAMPLE_FASTA = """\
>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT
"""


# ---------------------------------------------------------------------------
# HAMM -- Hamming distance
# ---------------------------------------------------------------------------

def hamming(a, b):
    """Mismatch count. Requires equal length and permits no gaps -- which is
    exactly the limitation that forces us into alignment (see alignment.py)."""
    if len(a) != len(b):
        raise ValueError("Hamming distance needs equal-length strings")
    return sum(x != y for x, y in zip(a, b))


# ---------------------------------------------------------------------------
# SUBS -- motif locations
# ---------------------------------------------------------------------------

def find_motif(seq, motif):
    """All 1-based start positions, INCLUDING overlapping matches.

    Overlaps matter biologically (repeats, tandem arrays) and are exactly what
    a naive regex findall would silently drop.
    """
    positions, start = [], 0
    while True:
        i = seq.find(motif, start)
        if i == -1:
            return positions
        positions.append(i + 1)   # Rosalind and genomics are 1-based
        start = i + 1


# ---------------------------------------------------------------------------
# CONS -- consensus and profile matrix
# ---------------------------------------------------------------------------

def profile_matrix(seqs):
    """Count each base at each column. This is a position weight matrix, the
    foundation of motif finding -- how you describe 'a transcription factor
    binds roughly TATAAA' as a probabilistic object rather than a fixed string.
    """
    n = len(seqs[0])
    profile = {b: [0] * n for b in "ACGT"}
    for seq in seqs:
        for i, base in enumerate(seq):
            if base in profile:
                profile[base][i] += 1
    return profile


def consensus(profile):
    n = len(next(iter(profile.values())))
    return "".join(max("ACGT", key=lambda b: profile[b][i]) for i in range(n))


# ---------------------------------------------------------------------------

def main():
    rule = "-" * 72
    records = list(parse_fasta(SAMPLE_FASTA))

    print(rule)
    print("GC  -- parse FASTA, find the highest-GC sequence")
    print(rule)
    for header, seq in records:
        print(f"  {header}  {len(seq):>3} bp   GC {gc_content(seq):.4%}")
    best = max(records, key=lambda r: gc_content(r[1]))
    print(f"  -> highest GC: {best[0]}")
    print("  GC content varies by organism and by region; it also predicts")
    print("  sequencing difficulty, which is why QC tools always report it.")

    print()
    print(rule)
    print("HAMM -- counting mismatches")
    print(rule)
    a, b = "GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"
    print(f"  {a}")
    print(f"  {''.join('|' if x == y else ' ' for x, y in zip(a, b))}")
    print(f"  {b}")
    print(f"  distance = {hamming(a, b)}")
    print("  Note what this CANNOT do: if b had one base deleted, every")
    print("  position after it would count as a mismatch. That failure is")
    print("  the entire motivation for alignment.")

    print()
    print(rule)
    print("SUBS -- motif occurrences, overlaps included")
    print(rule)
    seq, motif = "GATATATGCATATACTT", "ATAT"
    print(f"  sequence {seq}")
    print(f"  motif    {motif}")
    hits = find_motif(seq, motif)
    print(f"  1-based positions: {hits}")
    marks = [" "] * len(seq)
    for h in hits:
        marks[h - 1] = "^"
    print(f"           {''.join(marks)}")
    print("  positions 2 and 4 overlap -- a regex findall would miss one.")

    print()
    print(rule)
    print("PROT -- translation, reusing module 01")
    print(rule)
    rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
    print(f"  RNA     {rna}")
    print(f"  protein {translate(rna)}")

    print()
    print(rule)
    print("REVC -- and why it is everywhere")
    print(rule)
    s = "AAAACCCGGT"
    print(f"  {s}  ->  {reverse_complement(s)}")

    print()
    print(rule)
    print("CONS -- profile matrix and consensus")
    print(rule)
    seqs = ["ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC", "TTGGAACT",
            "ATGCCATT", "ATGGCACT"]
    for s in seqs:
        print(f"    {s}")
    prof = profile_matrix(seqs)
    print()
    for base in "ACGT":
        print(f"  {base}: {' '.join(f'{c}' for c in prof[base])}")
    print(f"  consensus: {consensus(prof)}")
    print("  Each column is a distribution, not a letter. Real motif models")
    print("  (TF binding sites, splice sites) are exactly this, log-scaled.")


if __name__ == "__main__":
    main()
