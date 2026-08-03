"""Footgun #1, made concrete: 0-based vs 1-based coordinates.

Read this once, believe it, and you will still get bitten -- but you will
recognise the symptom, which is the actual goal.

Run: python 03-formats/coordinates.py
"""

RULE = "-" * 74

SEQ = "ACGTACGTAC"          # a 10-base toy chromosome


def main():
    print(RULE)
    print("The same 4-base feature, written in every format")
    print(RULE)
    print(f"  sequence      {SEQ}")
    print("  index (0-based) 0123456789")
    print("  index (1-based) 1234567890")
    print()
    print("  We want to describe the GTAC at positions 3-6 (1-based).")
    print()
    print("  BED     chr1  2  6      0-based, half-open  -> start=2, end=6")
    print("  GFF     chr1  3  6      1-based, closed     -> start=3, end=6")
    print("  VCF     POS=3           1-based")
    print("  SAM     POS=3           1-based in the file")
    print("  pysam   start=2         0-based in the API, reading the same file")
    print()
    print(f"  Python slice:  SEQ[2:6] = {SEQ[2:6]!r}   <- BED numbers work directly")
    print(f"  Naive 1-based: SEQ[3:6] = {SEQ[3:6]!r}    <- off by one, silently")

    print()
    print(RULE)
    print("Why BED chose 0-based half-open (it is not arbitrary)")
    print(RULE)
    print("  length          = end - start        (no +1 to forget)")
    print("  adjacent blocks = [0,5) and [5,10)   (no ambiguity about 5)")
    print("  empty interval  = [5,5)              (representable at all)")
    print("  It is the same reason Python slices work that way. GFF predates")
    print("  the argument and stayed 1-based because biologists count from 1.")

    print()
    print(RULE)
    print("The conversion, and the two directions people get wrong")
    print(RULE)
    print("  BED -> 1-based:   start + 1, end unchanged")
    print("  1-based -> BED:   start - 1, end unchanged")
    print("  Only the START moves. Adjusting the end too is the classic bug,")
    print("  and it produces intervals one base short -- which looks fine")
    print("  until a variant sits on the boundary.")
    print()
    for bed_start, bed_end in [(2, 6), (0, 10), (5, 5)]:
        one_start, one_end = bed_start + 1, bed_end
        length = bed_end - bed_start
        print(f"  BED [{bed_start},{bed_end})  ==  1-based {one_start}-{one_end}"
              f"   length {length}   seq {SEQ[bed_start:bed_end]!r}")

    print()
    print(RULE)
    print("Footgun #2: the same position in two reference builds")
    print(RULE)
    print("  rs334, the sickle-cell variant:")
    print("    hg19 / GRCh37   chr11:5,248,232")
    print("    hg38 / GRCh38   chr11:5,227,002")
    print("  A 21,230-base difference. Neither is wrong; they are different")
    print("  coordinate systems for the same molecule, because hg38 fixed")
    print("  assembly errors and added sequence upstream.")
    print()
    print("  What makes this dangerous is that it does not error. Feed hg19")
    print("  coordinates to an hg38 reference and you get a real base at a")
    print("  real position in a real gene -- just the wrong one. The only")
    print("  reliable tripwire is asserting that the REF allele in your VCF")
    print("  actually matches the reference base (make_data.py does this).")
    print()
    print("  And the naming split, which is pure gratuitous pain:")
    print("    UCSC     chr11   chrM    chrX")
    print("    Ensembl  11      MT      X")
    print("  Tools do not error on a mismatch -- they find zero overlaps and")
    print("  report an empty result, which reads as a biological finding.")

    print()
    print(RULE)
    print("Habits that make these bugs cheap instead of expensive")
    print(RULE)
    print("  1. Put the build in every filename:   sample.hg38.bam")
    print("  2. Assert REF matches the reference before trusting a VCF.")
    print("  3. When a tool returns zero results, suspect naming before biology.")
    print("  4. Never hand-convert coordinates; use liftOver/CrossMap and")
    print("     accept that some regions have no equivalent at all.")
    print("  5. Look at the BAM in IGV. Off-by-one errors are visible there")
    print("     in a way they never are in a summary statistic.")


if __name__ == "__main__":
    main()
