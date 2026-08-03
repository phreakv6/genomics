"""The central dogma, implemented from scratch.

Every function here corresponds to a section of theory.md. No dependencies --
the point is that DNA->RNA->protein is string manipulation once you see it.

Run: python 01-central-dogma/central_dogma.py
"""

# ---------------------------------------------------------------------------
# 1. DNA as a string
# ---------------------------------------------------------------------------

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def count_bases(seq):
    """The literal 'hello world' of bioinformatics (Rosalind problem DNA)."""
    return {base: seq.count(base) for base in "ACGT"}


def gc_content(seq):
    """Fraction of G+C. Matters because GC-rich regions bind more tightly, are
    harder to sequence, and GC% is a species/region fingerprint."""
    return (seq.count("G") + seq.count("C")) / len(seq)


def complement(seq):
    return "".join(COMPLEMENT[b] for b in seq)


def reverse_complement(seq):
    """The other strand, written 5'->3' as convention demands.

    This is the most-used operation in the field: roughly half of all sequencer
    reads align to the reverse strand, so every aligner does this constantly.
    """
    return "".join(COMPLEMENT[b] for b in reversed(seq))


# ---------------------------------------------------------------------------
# 2. Transcription: DNA -> RNA
# ---------------------------------------------------------------------------

def transcribe(dna):
    """Thymine becomes uracil. That is the entire operation at string level."""
    return dna.replace("T", "U")


# ---------------------------------------------------------------------------
# 3. Translation: RNA -> protein
# ---------------------------------------------------------------------------

# The standard genetic code. 64 codons -> 20 amino acids + stop.
# Note how often the third base is irrelevant -- that degeneracy is why many
# DNA mutations are "synonymous" and change nothing about the protein.
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",   # AUG = Met = START
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",   # * = STOP
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Three-letter names, for reading variant notation like p.Glu6Val.
AA_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu",
    "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Ter",
}


def translate(rna, stop_at_stop=True):
    """Read three bases at a time through the codon table."""
    protein = []
    for i in range(0, len(rna) - len(rna) % 3, 3):
        aa = CODON_TABLE[rna[i:i + 3]]
        if aa == "*" and stop_at_stop:
            break
        protein.append(aa)
    return "".join(protein)


def six_frame_translation(dna):
    """Three frames forward, three on the reverse complement.

    When you don't know where a gene starts, this is how you look.
    """
    frames = {}
    for offset in range(3):
        frames[f"+{offset + 1}"] = translate(transcribe(dna[offset:]), stop_at_stop=False)
    rc = reverse_complement(dna)
    for offset in range(3):
        frames[f"-{offset + 1}"] = translate(transcribe(rc[offset:]), stop_at_stop=False)
    return frames


STOPS = {"UAA", "UAG", "UGA"}


def find_orfs(dna, min_aa=20):
    """Open reading frames: AUG ... stop, in any of the six frames.

    Coordinates are always reported against the ORIGINAL forward sequence,
    0-based half-open [start, end) like a BED interval, with `end` including
    the stop codon. This matters: a minus-strand ORF is found by searching the
    reverse complement, so its offset in that searched string counts from the
    far end of the input. Reporting that raw offset would make plus- and
    minus-strand hits incomparable, and would place the feature at the wrong
    end of the sequence. Mapping back is the `L - end` arithmetic below.

    Naive on purpose. Real gene finding has to deal with introns, which this
    cannot see -- that is exactly why annotation files exist (module 03).
    """
    L = len(dna)
    orfs = []
    for strand, strand_seq in (("+", dna), ("-", reverse_complement(dna))):
        rna = transcribe(strand_seq)
        for start in range(len(rna) - 2):
            if rna[start:start + 3] != "AUG":
                continue
            peptide = translate(rna[start:], stop_at_stop=True)
            # translate() stops before the stop codon, so the codon that
            # follows the peptide must itself be a stop for this to be a
            # complete ORF rather than one running off the end of the input.
            stop_at = start + 3 * len(peptide)
            if rna[stop_at:stop_at + 3] not in STOPS:
                continue
            if len(peptide) < min_aa:
                continue

            end = stop_at + 3                      # include the stop codon
            if strand == "+":
                fwd_start, fwd_end = start, end
            else:
                # index i of the reverse complement is index L-1-i of the
                # forward sequence, so the interval flips end-for-end.
                fwd_start, fwd_end = L - end, L - start

            orfs.append({
                "strand": strand,
                "start": fwd_start,                # 0-based, on the input seq
                "end": fwd_end,                    # half-open, includes stop
                "peptide": peptide,
            })
    return sorted(orfs, key=lambda o: (o["start"], o["strand"]))


def extract_orf(dna, orf):
    """Pull an ORF's nucleotides back out of the original sequence.

    Round-tripping through this is how you prove the coordinates are right --
    translating the extracted slice must reproduce the reported peptide.
    """
    chunk = dna[orf["start"]:orf["end"]]
    if orf["strand"] == "-":
        chunk = reverse_complement(chunk)
    return chunk


# ---------------------------------------------------------------------------
# 4. A real gene, and a real mutation
# ---------------------------------------------------------------------------

# First 31 codons of the human HBB gene (beta-globin) coding sequence,
# chromosome 11. Picked because the most famous point mutation in medicine
# lives at codon 6 of the mature protein.
HBB_CDS = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGAT"
    "GAAGTTGGTGGTGAGGCCCTGGGCAGG"
)

# A constructed sequence for demonstrating ORF finding on both strands.
# Built as:  filler + [forward ORF] + filler + [revcomp of a second ORF] + filler
# Real bacterial genomes genuinely look a bit like this -- genes packed on
# both strands with short gaps. Eukaryotic DNA does not, because of introns.
_ORF_FWD = "ATG" + "GCTTTAGGTCATAAGCCTTGGGAAGATTTC" + "TAA"   # 12 aa + stop
_ORF_REV = "ATG" + "AAACGTGTTTGGCAAGACTTAATTCCGAGT" + "TGA"   # 12 aa + stop
DEMO = ("CCTTACGAGT"
        + _ORF_FWD
        + "GGATCCTTAA"
        + "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[b]
                 for b in reversed(_ORF_REV))          # planted on minus strand
        + "TTAGGCCATA")


def apply_snv(seq, zero_based_pos, ref, alt):
    """Apply a single-nucleotide change, asserting the reference matches.

    That assertion is not paranoia -- 'my REF doesn't match the reference'
    is the classic symptom of a genome-build mismatch (hg19 vs hg38).
    """
    assert seq[zero_based_pos] == ref, (
        f"expected {ref} at {zero_based_pos}, found {seq[zero_based_pos]}"
    )
    return seq[:zero_based_pos] + alt + seq[zero_based_pos + 1:]


def classify_snv(cds, pos, ref, alt):
    """Sort a coding SNV into synonymous / missense / nonsense.

    This is the first and most basic step of variant interpretation.
    """
    mutant = apply_snv(cds, pos, ref, alt)
    wt_protein = translate(transcribe(cds), stop_at_stop=False)
    mut_protein = translate(transcribe(mutant), stop_at_stop=False)

    codon_index = pos // 3
    wt_aa, mut_aa = wt_protein[codon_index], mut_protein[codon_index]

    if wt_aa == mut_aa:
        kind = "synonymous"
    elif mut_aa == "*":
        kind = "nonsense"
    else:
        kind = "missense"

    # Protein numbering conventionally drops the initiator methionine,
    # which is why the sickle variant is called Glu6Val and not Glu7Val.
    protein_pos = codon_index  # 0-based codon index == 1-based mature position
    hgvs = f"p.{AA_NAMES[wt_aa]}{protein_pos}{AA_NAMES[mut_aa]}"
    return kind, hgvs, wt_protein, mut_protein


# ---------------------------------------------------------------------------

def main():
    seq = HBB_CDS
    rule = "-" * 72

    print(rule)
    print("1. DNA as a string -- human HBB (beta-globin), first 31 codons")
    print(rule)
    print(f"  sequence   {seq[:60]}...")
    print(f"  length     {len(seq)} bp  ({len(seq) // 3} codons)")
    print(f"  counts     {count_bases(seq)}")
    print(f"  GC content {gc_content(seq):.1%}")
    print()
    print(f"  forward    5'-{seq[:30]}...-3'")
    print(f"  complement 3'-{complement(seq)[:30]}...-5'")
    print(f"  rev-comp   5'-{reverse_complement(seq)[:30]}...-3'")
    print("  (rev-comp is the OTHER strand written the conventional way)")

    print()
    print(rule)
    print("2. Transcription: DNA -> RNA (T becomes U)")
    print(rule)
    rna = transcribe(seq)
    print(f"  DNA  {seq[:48]}")
    print(f"  RNA  {rna[:48]}")

    print()
    print(rule)
    print("3. Translation: RNA -> protein (three bases at a time)")
    print(rule)
    protein = translate(rna)
    for i in range(0, 24, 12):
        codons = " ".join(rna[j:j + 3] for j in range(i * 3, i * 3 + 36, 3))
        aas = "   ".join(protein[i:i + 12])
        print(f"  {codons}")
        print(f"    {aas}")
    print(f"\n  protein  {protein}")
    print("  ^ that is the real N-terminus of human haemoglobin beta chain")

    print()
    print(rule)
    print("4. Reading frames -- where you start changes everything")
    print(rule)
    for label, peptide in six_frame_translation(seq).items():
        offset = int(label[1]) - 1
        usable = len(seq) - offset
        leftover = usable % 3
        note = f"{len(peptide):>2} aa" + (f"  (+{leftover} base{'s' * (leftover > 1)} "
                                          "left over, discarded)" if leftover else "")
        print(f"  frame {label}  {peptide}  {note}")
    print()
    print("  only frame +1 is the real one; the rest are riddled with stops (*)")
    print()
    print("  Note the lengths differ. Frames +2/+3 start 1 or 2 bases in, so the")
    print("  sequence no longer divides evenly by 3 and the trailing 1-2 bases")
    print("  cannot form a codon. translate() drops them, via the")
    print("  'len(rna) - len(rna) % 3' bound on its loop. Two bases are not a")
    print("  codon and nothing can be said about them -- so with a 93bp input")
    print("  you get 31, 30, 30 residues, not 31 three times.")

    print()
    print(rule)
    print("5. A real mutation: sickle cell disease")
    print(rule)
    # HBB codon 6 (mature numbering): GAG -> GTG, a single A->T at CDS pos 20
    # (1-based). This one base is the entire molecular basis of sickle cell.
    kind, hgvs, wt, mut = classify_snv(seq, 19, "A", "T")
    print(f"  change     c.20A>T   (GAG -> GTG)")
    print(f"  effect     {kind}   {hgvs}")
    print(f"  wild type  {wt[:16]}")
    print(f"  mutant     {mut[:16]}")
    print("                   ^ Glu (charged, water-loving) -> Val (greasy)")
    print("  That hydrophobic patch makes haemoglobin molecules stick to each")
    print("  other and polymerise, deforming the red cell into a sickle.")
    print("  One base out of 3 billion.")

    print()
    print("  For contrast, a third-position change in the same codon:")
    kind2, hgvs2, _, _ = classify_snv(seq, 20, "G", "A")   # GAG -> GAA
    print(f"  change     c.21G>A   (GAG -> GAA)   effect: {kind2}  {hgvs2}")
    print("  Same codon, different position, zero consequence. That asymmetry")
    print("  is the degeneracy of the genetic code, and it is why variant")
    print("  annotation is the first filter in every clinical pipeline.")

    print()
    print(rule)
    print("6. ORF finding (naive)")
    print(rule)
    orfs = find_orfs(seq, min_aa=5)
    if orfs:
        for orf in orfs:
            print(f"  {orf['strand']} {orf['start']:>3}-{orf['end']:<3} {orf['peptide']}")
    else:
        print("  Nothing found in HBB_CDS: this fragment is cut mid-gene, so the")
        print("  real ORF runs off the end without ever reaching a stop codon.")
        print("  Realistic -- genes do not come pre-trimmed, and introns make")
        print("  naive ORF finding fail outright on real genomic DNA.")

    print()
    print("  So here is a constructed sequence with one ORF on each strand:")
    print()
    orfs = find_orfs(DEMO, min_aa=5)
    print(f"  sequence ({len(DEMO)} bp)  {DEMO[:56]}...")
    print()
    print(f"  {'strand':<8}{'start':>6}{'end':>6}{'len':>6}   peptide")
    for orf in orfs:
        nt = orf["end"] - orf["start"]
        print(f"  {orf['strand']:<8}{orf['start']:>6}{orf['end']:>6}{nt:>5}nt   "
              f"{orf['peptide']}*")

    print()
    print("  Coordinates are on the FORWARD sequence, 0-based half-open, with")
    print("  the stop codon included -- so they can be compared to each other,")
    print("  and dropped straight into a BED file.")
    print()
    print("  Proof they are right: slice the original sequence at those")
    print("  coordinates, reverse-complement if the hit was on the minus")
    print("  strand, and translate. It must reproduce the peptide.")
    print()
    for orf in orfs:
        nts = extract_orf(DEMO, orf)
        back = translate(transcribe(nts), stop_at_stop=True)
        ok = "OK" if back == orf["peptide"] else "MISMATCH"
        print(f"    strand {orf['strand']}  DEMO[{orf['start']}:{orf['end']}]"
              f"{' -> revcomp' if orf['strand'] == '-' else '           '}"
              f" -> {back}   [{ok}]")
        assert back == orf["peptide"], "coordinate mapping is wrong"
    print()
    print("  Note the minus-strand ORF's start is NOT its offset in the")
    print("  reverse-complemented string that was searched -- that offset")
    print("  counts from the opposite end. Mapping back is 'L - end', and")
    print("  getting it wrong silently places the feature at the wrong end")
    print("  of the sequence, which is the reading-frame version of the")
    print("  coordinate bugs in module 03.")


if __name__ == "__main__":
    main()
