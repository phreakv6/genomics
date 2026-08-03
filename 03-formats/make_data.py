"""Fetch a real slice of the human genome and simulate reads from it.

Region: chr11:5,225,464-5,229,395 (hg38) -- the human HBB gene, the same gene
whose first 31 codons module 01 translated by hand.

We simulate reads rather than download real ones so that the ground truth is
known exactly: we plant the sickle-cell variant (rs334, chr11:5227002 T>A) as a
heterozygote, then see whether the pipeline in run.sh finds it. Module 04 does
the same thing with real reads and a published truth set.

Run: python 03-formats/make_data.py
"""

import json
import os
import random
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# hg38 coordinates. Note the build in the variable name -- adopt this habit.
CHROM = "chr11"
START_1BASED = 5225464          # 1-based inclusive, as a human would quote it
END_1BASED = 5229395
BUILD = "hg38"

# rs334: the sickle-cell variant. On the PLUS strand of hg38 it is T>A.
# HBB is transcribed from the MINUS strand, which is why module 01 saw it as
# an A>T change in the coding sequence. Same event, two reading conventions.
RS334_POS = 5227002             # 1-based
RS334_REF, RS334_ALT = "T", "A"

COVERAGE = 50
READ_LEN = 150
ERROR_RATE = 0.001              # ~Q30
random.seed(11)                 # reproducible: pipelines must be re-runnable


def fetch_reference():
    """UCSC's REST API. Note it takes 0-BASED, half-open coordinates, while
    the numbers above are 1-based -- so we subtract 1 from the start only.
    This conversion is footgun #1 in miniature."""
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome={BUILD}"
           f";chrom={CHROM};start={START_1BASED - 1};end={END_1BASED}")
    with urllib.request.urlopen(url) as r:
        return json.load(r)["dna"].upper()


def write_fasta(path, name, seq, width=60):
    """FASTA: a '>' header, then sequence wrapped. The wrapping is why you
    cannot treat one line as one record."""
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def phred_char(q):
    """Phred+33 encoding: Q = -10*log10(P(error)), stored as chr(Q+33)."""
    return chr(min(q, 40) + 33)


def simulate_reads(ref, ref_offset, fastq_path):
    """Two haplotypes: one reference, one carrying rs334. Sample reads from
    both, add sequencing errors, write FASTQ."""
    hap_ref = ref
    idx = RS334_POS - ref_offset            # 0-based index into our slice
    assert hap_ref[idx] == RS334_REF, "reference base mismatch -- wrong build?"
    hap_alt = hap_ref[:idx] + RS334_ALT + hap_ref[idx + 1:]

    n_reads = (len(ref) * COVERAGE) // READ_LEN
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

    with open(fastq_path, "w") as f:
        for i in range(n_reads):
            hap = hap_ref if i % 2 == 0 else hap_alt     # heterozygous: 50/50
            start = random.randint(0, len(hap) - READ_LEN)
            read = list(hap[start:start + READ_LEN])
            quals = []
            for j, base in enumerate(read):
                if random.random() < ERROR_RATE:
                    read[j] = random.choice([b for b in "ACGT" if b != base])
                    quals.append(random.randint(15, 25))  # errors skew low-Q
                else:
                    quals.append(random.randint(30, 40))
            seq = "".join(read)
            # Half of all real reads come off the reverse strand.
            if i % 4 in (2, 3):
                seq = "".join(complement[b] for b in reversed(seq))
                quals = quals[::-1]
            qual = "".join(phred_char(q) for q in quals)
            # The 4 lines of FASTQ, always in this order:
            f.write(f"@read{i}\n{seq}\n+\n{qual}\n")

    return n_reads


def main():
    os.makedirs(DATA, exist_ok=True)
    ref_path = os.path.join(DATA, f"HBB.{BUILD}.fa")
    fastq_path = os.path.join(DATA, f"sample.{BUILD}.fastq")

    print(f"fetching {CHROM}:{START_1BASED:,}-{END_1BASED:,} ({BUILD}) from UCSC...")
    ref = fetch_reference()
    print(f"  got {len(ref):,} bp")

    # Name the contig for the slice, and record the offset so coordinates can
    # be translated back to real chromosome positions later.
    contig = f"{CHROM}:{START_1BASED}-{END_1BASED}"
    write_fasta(ref_path, contig, ref)
    print(f"  wrote {ref_path}")

    n = simulate_reads(ref, START_1BASED, fastq_path)
    print(f"simulated {n:,} reads of {READ_LEN}bp at ~{COVERAGE}x coverage")
    print(f"  planted: {CHROM}:{RS334_POS} {RS334_REF}>{RS334_ALT} "
          f"heterozygous (rs334, sickle cell)")
    print(f"  wrote {fastq_path}")
    print()
    print("The planted variant is the ground truth. run.sh must recover it.")


if __name__ == "__main__":
    main()
