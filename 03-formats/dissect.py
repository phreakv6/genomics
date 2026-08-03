"""Take apart the files run.sh produced, field by field.

Run: python 03-formats/dissect.py   (after run.sh)
"""

import os
import pysam

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
RULE = "-" * 74
CONTIG_START = 5225464          # our slice's 1-based start on chr11


def phred_to_error(q):
    return 10 ** (-q / 10)


def fastq_quality():
    print(RULE); print("FASTQ -- the quality string is the whole point"); print(RULE)
    with open(os.path.join(DATA, "sample.hg38.fastq")) as f:
        name, seq, _, qual = [next(f).rstrip() for _ in range(4)]
    print(f"  1  {name}          <- identifier")
    print(f"  2  {seq[:50]}...   <- bases called")
    print(f"  3  +                    <- separator (historically a repeated id)")
    print(f"  4  {qual[:50]}...   <- one quality char per base")
    print()
    print("  decoding the first few:  Q = ord(char) - 33")
    print(f"    {'char':>6} {'Q':>4} {'P(error)':>12}   meaning")
    for ch in qual[:4] + "".join(sorted(set(qual))[:2]):
        q = ord(ch) - 33
        note = "1 in " + f"{1/phred_to_error(q):,.0f}"
        print(f"    {ch:>6} {q:>4} {phred_to_error(q):>12.6f}   {note}")
    print()
    print("  Q30 = 1 error in 1,000 = the usual 'good data' threshold.")
    print("  These numbers are not metadata -- the variant caller in step 5")
    print("  multiplies them together to decide if a pile of bases is a real")
    print("  variant or noise. Bad quality does not just get filtered; it gets")
    print("  weighted. This is why variant calling is Bayesian.")


def decode_flag(flag):
    bits = [
        (0x1, "paired"), (0x2, "properly paired"), (0x4, "unmapped"),
        (0x8, "mate unmapped"), (0x10, "REVERSE strand"), (0x20, "mate reverse"),
        (0x40, "first in pair"), (0x80, "second in pair"),
        (0x100, "secondary alignment"), (0x200, "QC fail"),
        (0x400, "PCR duplicate"), (0x800, "supplementary"),
    ]
    return [name for bit, name in bits if flag & bit]


def sam_fields():
    print(); print(RULE); print("SAM/BAM -- 11 columns, two of them encoded"); print(RULE)
    bam = pysam.AlignmentFile(os.path.join(DATA, "sorted.bam"), "rb")
    reads = []
    for r in bam:
        reads.append(r)
        if len(reads) >= 400:
            break
    fwd = next(r for r in reads if not r.is_reverse)
    rev = next(r for r in reads if r.is_reverse)

    for label, r in (("forward-strand read", fwd), ("reverse-strand read", rev)):
        print(f"\n  {label}: {r.query_name}")
        print(f"    FLAG   {r.flag:>6}  = {' + '.join(decode_flag(r.flag)) or 'all bits clear'}")
        print(f"    RNAME  {r.reference_name}")
        # pysam is 0-based; the SAM file on disk is 1-based. This is footgun #1
        # hiding inside a library API.
        print(f"    POS    {r.reference_start + 1:>6}  (pysam gave {r.reference_start}, "
              f"0-based -- the file says {r.reference_start + 1})")
        print(f"    MAPQ   {r.mapping_quality:>6}  (60 = uniquely placed; 0 = maps elsewhere too)")
        print(f"    CIGAR  {r.cigarstring:>6}  = {explain_cigar(r.cigarstring)}")

    print()
    print("  MAPQ distribution across the first 400 reads:")
    from collections import Counter
    for mq, n in sorted(Counter(r.mapping_quality for r in reads).items(), reverse=True):
        print(f"    MAPQ {mq:>3}  {'#' * (n * 40 // len(reads)):<40} {n}")
    print("  All 60 here because this is one clean, non-repetitive gene.")
    print("  On a whole genome you would see a fat MAPQ 0 spike -- repeats,")
    print("  where a read genuinely cannot be placed. Filtering those out is")
    print("  the highest-leverage single filter in the field.")
    bam.close()


def explain_cigar(cigar):
    import re
    names = {"M": "aligned", "I": "insertion", "D": "deletion", "N": "skipped",
             "S": "soft-clipped", "H": "hard-clipped", "P": "padding",
             "=": "match", "X": "mismatch"}
    parts = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    return ", ".join(f"{n} {names[op]}" for n, op in parts)


def vcf_fields():
    print(); print(RULE); print("VCF -- a diff against the reference, not your DNA"); print(RULE)
    vcf = pysam.VariantFile(os.path.join(DATA, "calls.vcf"))
    for rec in vcf:
        real_pos = CONTIG_START + rec.pos - 1
        print(f"\n  CHROM   {rec.chrom}")
        print(f"  POS     {rec.pos}   (contig-local; real coordinate chr11:{real_pos:,})")
        print(f"  REF     {rec.ref}      <- what the reference genome says")
        print(f"  ALT     {','.join(str(a) for a in rec.alts)}      <- what this sample has")
        print(f"  QUAL    {rec.qual:.1f}  (Phred: P(no variant here) = "
              f"{phred_to_error(rec.qual):.2e})")
        dp = rec.info.get("DP")
        print(f"  INFO    DP={dp}  <- {dp} reads covered this position")
        for sample, values in rec.samples.items():
            gt = values["GT"]
            ad = values.get("AD")
            print(f"  SAMPLE  {sample}")
            print(f"    GT    {gt[0]}/{gt[1]}   <- 0=REF allele, 1=ALT allele. "
                  f"{'heterozygous' if gt[0] != gt[1] else 'homozygous'}")
            if ad:
                print(f"    AD    {ad[0]},{ad[1]}  <- {ad[0]} reads showed {rec.ref}, "
                      f"{ad[1]} showed {rec.alts[0]}")
                print(f"          ratio {ad[1] / sum(ad):.0%} alt -- near 50%, "
                      "as a heterozygote should be")
    print()
    print(f"  This is rs334. Reference says T at chr11:{CONTIG_START + rec.pos - 1:,};")
    print("  this sample carries one T and one A. HBB is on the minus strand,")
    print("  so on the coding strand that is the A>T that module 01 turned into")
    print("  p.Glu6Val -- sickle cell. Same event, three coordinate conventions.")


def pileup_view():
    print(); print(RULE); print("The pileup -- what the caller actually sees"); print(RULE)
    bam = pysam.AlignmentFile(os.path.join(DATA, "sorted.bam"), "rb")
    ref = pysam.FastaFile(os.path.join(DATA, "HBB.hg38.fa"))
    contig = bam.references[0]
    target0 = 5227002 - CONTIG_START          # 0-based index of rs334

    for col in bam.pileup(contig, target0, target0 + 1, truncate=True,
                          min_base_quality=0):
        bases = [p.alignment.query_sequence[p.query_position]
                 for p in col.pileups if p.query_position is not None]
        from collections import Counter
        counts = Counter(bases)
        refbase = ref.fetch(contig, col.reference_pos, col.reference_pos + 1)
        print(f"  chr11:{CONTIG_START + col.reference_pos:,}  reference base {refbase}")
        print(f"  {col.nsegments} reads stacked here:")
        for base, n in counts.most_common():
            tag = "(reference)" if base == refbase else "(alternate)"
            print(f"    {base}  {'#' * n:<40} {n:>3}  {tag}")
        print()
        print("  Two bases, roughly 50/50, both at high quality. That is the")
        print("  signature of a heterozygous site. One base at 50/50 with LOW")
        print("  qualities would be noise; a 5% minority would be contamination")
        print("  or somatic mosaicism. The caller distinguishes these purely")
        print("  from the Phred scores in the FASTQ.")
    bam.close()


if __name__ == "__main__":
    fastq_quality()
    sam_fields()
    vcf_fields()
    pileup_view()
