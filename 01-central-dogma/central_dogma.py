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


def find_orfs(dna, min_aa=20):
    """Open reading frames: AUG ... stop, in any of the six frames.

    Naive on purpose. Real gene finding has to deal with introns, which this
    cannot see -- that is exactly why annotation files exist (module 03).
    """
    orfs = []
    for label, strand_seq in (("+", dna), ("-", reverse_complement(dna))):
        rna = transcribe(strand_seq)
        for start in range(len(rna) - 2):
            if rna[start:start + 3] != "AUG":
                continue
            peptide = translate(rna[start:], stop_at_stop=True)
            # translate() stopped early only if a stop codon was actually hit
            hit_stop = start + 3 * len(peptide) + 3 <= len(rna)
            if hit_stop and len(peptide) >= min_aa:
                orfs.append((label, start, peptide))
    return orfs


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
        print(f"  frame {label}  {peptide[:30]}")
    print("  only frame +1 is the real one; the rest are riddled with stops (*)")

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
        for strand, start, peptide in orfs:
            print(f"  strand {strand}  start {start:>3}  {len(peptide):>3} aa  {peptide[:40]}")
    else:
        print("  none found with a stop codon -- this fragment is cut mid-gene,")
        print("  so the real ORF runs off the end. Realistic: genes rarely come")
        print("  pre-trimmed, and introns make naive ORF finding fail on real DNA.")


if __name__ == "__main__":
    main()
