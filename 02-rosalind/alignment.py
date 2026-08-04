"""Sequence alignment by dynamic programming -- the algorithmic core.

Needleman-Wunsch (global) and Smith-Waterman (local), with traceback, on
sequences small enough that the matrices are printable and readable.

This is what `bwa` runs inside the small candidate windows its BWT index
identifies. Understanding this is understanding module 04.

Run: python 02-rosalind/alignment.py
"""

MATCH, MISMATCH, GAP = 1, -1, -2


def score(a, b):
    return MATCH if a == b else MISMATCH


def needleman_wunsch(a, b):
    """Global alignment: force both sequences to align end to end.

    F(i,j) = best score aligning a[:i] against b[:j].
    """
    n, m = len(a), len(b)
    F = [[0] * (m + 1) for _ in range(n + 1)]

    # First row/column: aligning against nothing means all gaps.
    for i in range(1, n + 1):
        F[i][0] = i * GAP
    for j in range(1, m + 1):
        F[0][j] = j * GAP

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            F[i][j] = max(
                F[i - 1][j - 1] + score(a[i - 1], b[j - 1]),   # diagonal: align
                F[i - 1][j] + GAP,                             # up:   gap in b
                F[i][j - 1] + GAP,                             # left: gap in a
            )

    # Traceback from the bottom-right corner -- global means we must end there.
    ai, bi, i, j = [], [], n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + score(a[i - 1], b[j - 1]):
            ai.append(a[i - 1]); bi.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and F[i][j] == F[i - 1][j] + GAP:
            ai.append(a[i - 1]); bi.append("-"); i -= 1
        else:
            ai.append("-"); bi.append(b[j - 1]); j -= 1

    return F[n][m], "".join(reversed(ai)), "".join(reversed(bi)), F


def smith_waterman(a, b):
    """Local alignment: find the best-matching subregion.

    Two changes from Needleman-Wunsch:
      1. clamp negative scores to 0 (a bad prefix is abandoned, not carried)
      2. start traceback from the maximum cell anywhere, stop when you hit 0
    """
    n, m = len(a), len(b)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    best, best_pos = 0, (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            F[i][j] = max(
                0,                                             # <- the change
                F[i - 1][j - 1] + score(a[i - 1], b[j - 1]),
                F[i - 1][j] + GAP,
                F[i][j - 1] + GAP,
            )
            if F[i][j] > best:
                best, best_pos = F[i][j], (i, j)

    ai, bi = [], []
    i, j = best_pos
    while i > 0 and j > 0 and F[i][j] > 0:
        if F[i][j] == F[i - 1][j - 1] + score(a[i - 1], b[j - 1]):
            ai.append(a[i - 1]); bi.append(b[j - 1]); i -= 1; j -= 1
        elif F[i][j] == F[i - 1][j] + GAP:
            ai.append(a[i - 1]); bi.append("-"); i -= 1
        else:
            ai.append("-"); bi.append(b[j - 1]); j -= 1

    return best, "".join(reversed(ai)), "".join(reversed(bi)), F, best_pos


def print_matrix(F, a, b, highlight=None):
    print("        " + "".join(f"{c:>5}" for c in "-" + b))
    for i, row in enumerate(F):
        label = "-" if i == 0 else a[i - 1]
        cells = []
        for j, v in enumerate(row):
            s = f"{v:>5}"
            if highlight == (i, j):
                s = f"{('[' + str(v) + ']'):>5}"
            cells.append(s)
        print(f"   {label:>2}   " + "".join(cells))


def show(a, b, ali_a, ali_b):
    bars = "".join("|" if x == y else (" " if "-" in (x, y) else ".")
                   for x, y in zip(ali_a, ali_b))
    print(f"    {ali_a}")
    print(f"    {bars}   | match   . mismatch   (blank) gap")
    print(f"    {ali_b}")
    matches = bars.count("|")
    print(f"    {matches}/{len(ali_a)} identical  "
          f"({bars.count('.')} mismatches, {ali_a.count('-') + ali_b.count('-')} gaps)")


def main():
    rule = "-" * 72

    a, b = "GATTACA", "GCATGCU"

    print(rule)
    print("NEEDLEMAN-WUNSCH -- global alignment")
    print(rule)
    print(f"  scoring: match {MATCH:+d}, mismatch {MISMATCH:+d}, gap {GAP:+d}")
    print(f"  a = {a}   b = {b}\n")
    s, aa, bb, F = needleman_wunsch(a, b)
    print_matrix(F, a, b)
    print(f"\n  score {s}, traceback from the bottom-right corner:\n")
    show(a, b, aa, bb)
    print("\n  Note the first row and column: aligning against nothing costs")
    print("  one gap penalty per character. Global alignment has no choice but")
    print("  to consume both sequences entirely.")

    print()
    print(rule)
    print("SMITH-WATERMAN -- local alignment")
    print(rule)
    # A short 'read' that matches the middle of a longer 'reference', with one
    # substitution -- the situation every sequencer read is in.
    ref = "TTACGGATCAGCTTAGCATCGGCTAAGTTCAG"
    read = "GATCAGCTTAGCTTCGG"
    print(f"  reference {ref}")
    print(f"  read      {read}   (matches the middle, with one mismatch)\n")
    s2, aa2, bb2, F2, pos = smith_waterman(ref, read)
    print(f"  best score {s2} at cell {pos}\n")
    show(ref, read, aa2, bb2)
    print("\n  The read was placed correctly despite the mismatch, and without")
    print("  penalising the unmatched flanks of the reference. That is exactly")
    print("  what read alignment needs, and why SW is the right shape for it.")

    print()
    print(rule)
    print("WHY THIS DOESN'T SCALE")
    print(rule)
    human_genome = 3_100_000_000
    read_len, n_reads = 150, 500_000_000
    cells_one = human_genome * read_len
    print(f"  one {read_len}bp read vs the human genome:")
    print(f"    {human_genome:,} x {read_len} = {cells_one:,.0f} DP cells")
    print(f"  a typical 30x whole genome is ~{n_reads:,} reads:")
    print(f"    {cells_one * n_reads:.2e} cells")
    print("  At a generous 1e9 cells/sec that is ~7 billion years.")
    print()
    print("  The fix: index the genome once with a Burrows-Wheeler transform")
    print("  (the bzip2 transform, repurposed as a search index). Exact-match")
    print("  lookup then costs O(read length), independent of genome size --")
    print("  find candidate positions in microseconds, and run the DP above")
    print("  only inside those few small windows.")
    print()
    print("  That is bwa. The biggest speedup in genomics was a data-structures")
    print("  trick, not a biological insight.")


if __name__ == "__main__":
    main()
