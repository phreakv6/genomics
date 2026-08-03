"""Score our variant calls against the GIAB truth set, and look at every error.

The numbers matter less than the error inspection at the bottom: three false
positives you understand are worth more than a precision figure you don't.

Run: python 04-pipeline/benchmark.py   (after run.sh)
"""

import os
import pysam

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
RULE = "-" * 74


def load(path):
    """Return {(chrom, pos, ref, alt): record}, split by variant class."""
    out = {}
    for rec in pysam.VariantFile(path):
        for alt in rec.alts or ():
            out[(rec.chrom, rec.pos, rec.ref, alt)] = rec
    return out


def kind(key):
    _, _, ref, alt = key
    return "SNP" if len(ref) == 1 and len(alt) == 1 else "INDEL"


def metrics(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    return precision, recall, f1


def main():
    truth = load(os.path.join(DATA, "truth.hc.vcf.gz"))
    calls = load(os.path.join(DATA, "calls.hc.vcf.gz"))

    tp_keys = truth.keys() & calls.keys()
    fp_keys = calls.keys() - truth.keys()
    fn_keys = truth.keys() - calls.keys()

    print(RULE)
    print("HG002 (GIAB), chr11:5,205,000-5,295,000, ~31x, high-confidence only")
    print(RULE)
    print(f"  truth variants   {len(truth)}")
    print(f"  our calls        {len(calls)}")
    print()
    print(f"  {'class':<8}{'TP':>6}{'FP':>6}{'FN':>6}   {'precision':>10}{'recall':>9}{'F1':>8}")
    for cls in ("SNP", "INDEL", "ALL"):
        sel = (lambda k: True) if cls == "ALL" else (lambda k, c=cls: kind(k) == c)
        tp = sum(1 for k in tp_keys if sel(k))
        fp = sum(1 for k in fp_keys if sel(k))
        fn = sum(1 for k in fn_keys if sel(k))
        p, r, f = metrics(tp, fp, fn)
        print(f"  {cls:<8}{tp:>6}{fp:>6}{fn:>6}   {p:>10.4f}{r:>9.4f}{f:>8.4f}")

    print()
    print("  Reference points for a production short-read pipeline:")
    print("    SNPs    ~0.995 precision / ~0.995 recall")
    print("    Indels  ~0.99 precision / ~0.95 recall, far worse in repeats")
    print("  We used bcftools (a pileup caller) with no BQSR, no haplotype")
    print("  reassembly and no filtering, so landing below that is expected.")

    print()
    print(RULE)
    print("FALSE POSITIVES -- we called it, GIAB says it isn't there")
    print(RULE)
    if not fp_keys:
        print("  none")
    for key in sorted(fp_keys, key=lambda k: k[1]):
        rec = calls[key]
        s = rec.samples[0]
        ad = s.get("AD")
        gt = "/".join(str(a) for a in s["GT"])
        print(f"\n  chr11:{key[1]:,}  {key[2]}>{key[3]}  ({kind(key)})")
        print(f"    QUAL {rec.qual:.1f}   GT {gt}   DP {rec.info.get('DP')}"
              + (f"   AD {','.join(map(str, ad))}" if ad else ""))
        if ad and sum(ad) and ad[1] / sum(ad) < 0.25:
            print("    -> alt fraction well under 25%: too low for a real het.")
            print("       Classic noise/misalignment artefact; a depth+AF filter")
            print("       would remove it at almost no cost to recall.")
        elif rec.qual < 50:
            print("    -> low QUAL: a quality threshold would drop this.")
        else:
            print("    -> looks convincing on paper. Candidates: misalignment in")
            print("       a repeat, or a representation difference GIAB writes")
            print("       another way. Look at it in IGV before believing either.")

    print()
    print(RULE)
    print("FALSE NEGATIVES -- GIAB says it's there, we missed it")
    print(RULE)
    if not fn_keys:
        print("  none")
    for key in sorted(fn_keys, key=lambda k: k[1]):
        rec = truth[key]
        s = rec.samples[0]
        gt = "/".join(str(a) for a in s["GT"])
        print(f"\n  chr11:{key[1]:,}  {key[2]}>{key[3]}  ({kind(key)})   truth GT {gt}")
        # Did we have any call nearby? Near-misses are representation issues,
        # not true misses.
        near = [k for k in calls if abs(k[1] - key[1]) <= 10]
        if near:
            print("    -> we DID call something within 10bp:")
            for k in sorted(near, key=lambda k: k[1]):
                print(f"         chr11:{k[1]:,} {k[2]}>{k[3]}")
            print("       Same event, written differently. This is exactly the")
            print("       representation problem hap.py/vcfeval exist to solve;")
            print("       bcftools isec cannot see that these are the same.")
        else:
            print("    -> nothing called nearby: a genuine miss. Check depth and")
            print("       mapping quality at this position in IGV.")

    print()
    print(RULE)
    print("Coverage of the missed and spurious sites")
    print(RULE)
    bam = pysam.AlignmentFile(os.path.join(DATA, "hg002.markdup.bam"), "rb")
    for label, keys in (("FN", fn_keys), ("FP", fp_keys)):
        for key in sorted(keys, key=lambda k: k[1]):
            pos0 = key[1] - 1
            depths, mapqs = 0, []
            for col in bam.pileup("chr11", pos0, pos0 + 1, truncate=True,
                                  min_base_quality=0):
                depths = col.nsegments
                mapqs = [p.alignment.mapping_quality for p in col.pileups]
            mq0 = sum(1 for m in mapqs if m == 0)
            print(f"  {label}  chr11:{key[1]:,}   depth {depths:>3}   "
                  f"MAPQ0 reads {mq0}/{len(mapqs)}")
    bam.close()
    print()
    print("  Low depth explains misses; a pile of MAPQ 0 reads explains")
    print("  spurious calls -- those are reads that could have come from")
    print("  somewhere else in the genome, and often did.")

    print()
    print(RULE)
    print("Is the sickle-cell site in here?")
    print(RULE)
    rs334 = [k for k in truth if k[1] == 5227002] + [k for k in calls if k[1] == 5227002]
    if rs334:
        print(f"  called/known variant at chr11:5,227,002: {rs334}")
    else:
        print("  No variant at chr11:5,227,002 in either truth or our calls.")
        print("  HG002 carries the reference T on both chromosomes -- he does")
        print("  not have the sickle-cell allele, so there is nothing to call.")
        print("  Absence from a VCF means 'matches the reference', not 'no data'.")
        print("  That distinction is the most common misreading of a VCF.")


if __name__ == "__main__":
    main()
