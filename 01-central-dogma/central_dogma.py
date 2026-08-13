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

# Three-letter names, for reading variant notation like p.Glu7Val.
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
#
# Two DIFFERENT genes at two DIFFERENT positions, one on each strand:
#
#   0-10   10 ---------- 46   46-56   56 ---------- 92   92-102
#   filler  [ _ORF_FWD ]      filler  [ _ORF_REV ]       filler
#            + strand                  - strand
#
# They are not reverse complements of each other -- they are separate genes
# that happen to sit near each other. Real bacterial genomes look like this:
# genes densely packed on both strands with short gaps between. (A pair that
# WERE reverse complements would be an overlapping/ambisense gene, where one
# stretch of DNA encodes two proteins read in opposite directions. Real, but
# largely confined to viruses under extreme genome-size pressure, because
# every base then has to satisfy two coding constraints at once.)
#
# Eukaryotic DNA does not look like this either, because introns interrupt
# the coding sequence -- which is why naive ORF finding fails on human DNA.
#
# Each is 12 codons: ATG + 10 more + a stop, giving an 11-residue peptide.
_ORF_FWD = "ATG" + "GCTTTAGGTCATAAGCCTTGGGAAGATTTC" + "TAA"
_ORF_REV = "ATG" + "AAACGTGTTTGGCAAGACTTAATTCCGAGT" + "TGA"

# _ORF_REV is written as the gene READS -- 5'->3' along the minus strand. To
# plant it we need what the PLUS strand must say at that location, which is
# its reverse complement. This is the same conversion find_orfs() has to undo.
DEMO = ("CCTTACGAGT"
        + _ORF_FWD
        + "GGATCCTTAA"
        + reverse_complement(_ORF_REV)
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


def classify_snv(cds, pos, ref, alt, cds_offset=0):
    """Sort a coding SNV into synonymous / missense / nonsense / etc.

    `cds_offset` is where the coding sequence starts within `cds`. It defaults
    to 0 because HBB_CDS is already trimmed to the start codon -- but the
    parameter exists to make the dependency explicit. Every result below is
    relative to that choice; see demo_frame_dependence() for what happens when
    it is wrong.

    Returns a dict, because a variant consequence is genuinely several facts
    and squashing them into a string loses the ones that matter.
    """
    mutant = apply_snv(cds, pos, ref, alt)
    wt_protein = translate(transcribe(cds[cds_offset:]), stop_at_stop=False)
    mut_protein = translate(transcribe(mutant[cds_offset:]), stop_at_stop=False)

    codon_index = (pos - cds_offset) // 3            # 0-based
    wt_aa, mut_aa = wt_protein[codon_index], mut_protein[codon_index]

    if wt_aa == mut_aa:
        kind = "synonymous"
    elif codon_index == 0:
        # Destroying the start codon means translation cannot initiate here at
        # all. The protein is not "one residue different" -- it is absent, or
        # starts at some downstream AUG in a different frame.
        kind = "start-lost"
    elif wt_aa == "*":
        kind = "stop-lost"                           # translation runs on
    elif mut_aa == "*":
        kind = "nonsense"                            # premature stop
    else:
        kind = "missense"

    # HGVS protein numbering counts from the initiator methionine, which is
    # residue 1. So codon_index 6 (0-based) is p.7.
    hgvs_pos = codon_index + 1
    hgvs_p = f"p.{AA_NAMES[wt_aa]}{hgvs_pos}{AA_NAMES[mut_aa]}"

    # Much older literature uses MATURE numbering, which counts after the
    # initiator methionine has been cleaved off the finished protein -- so it
    # is one lower, and undefined for the initiator itself. This is why the
    # sickle variant is universally known as Glu6Val while ClinVar records it
    # as p.Glu7Val. Same variant, two conventions, endless confusion.
    mature_p = (f"p.{AA_NAMES[wt_aa]}{codon_index}{AA_NAMES[mut_aa]}"
                if codon_index >= 1 else None)

    return {
        "kind": kind,
        "hgvs_c": f"c.{pos - cds_offset + 1}{ref}>{alt}",   # 1-based on the CDS
        "hgvs_p": hgvs_p,
        "legacy_p": mature_p,
        "wt_protein": wt_protein,
        "mut_protein": mut_protein,
    }


def demo_frame_dependence(cds, pos, alt):
    """The same base change, read in each of the three possible frames.

    Answers the question 'is this really the same codon?' -- it is not. Codon
    boundaries are a property of the ANNOTATION, not of the base.
    """
    rows = []
    for offset in (0, 1, 2):
        codon_index = (pos - offset) // 3
        codon_start = offset + 3 * codon_index
        wt_codon = cds[codon_start:codon_start + 3]
        if len(wt_codon) < 3:
            continue
        i = pos - codon_start
        mut_codon = wt_codon[:i] + alt + wt_codon[i + 1:]
        wt_aa = CODON_TABLE[transcribe(wt_codon)]
        mut_aa = CODON_TABLE[transcribe(mut_codon)]
        if wt_aa == mut_aa:
            kind = "synonymous"
        elif wt_aa == "*":
            kind = "stop-lost"
        elif mut_aa == "*":
            kind = "nonsense"
        else:
            kind = "missense"
        rows.append({
            "offset": offset, "codon_pos": i + 1,
            "wt_codon": wt_codon, "mut_codon": mut_codon,
            "wt_aa": wt_aa, "mut_aa": mut_aa, "kind": kind,
        })
    return rows


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
    # A single A->T at CDS position 20 (1-based). This one base is the entire
    # molecular basis of sickle cell disease.
    r = classify_snv(seq, 19, "A", "T")
    print(f"  change     {r['hgvs_c']}   (GAG -> GTG)")
    print(f"  effect     {r['kind']}   {r['hgvs_p']}")
    print(f"  wild type  {r['wt_protein'][:16]}")
    print(f"  mutant     {r['mut_protein'][:16]}")
    print("                   ^ Glu (charged, water-loving) -> Val (greasy)")
    print("  That hydrophobic patch makes haemoglobin molecules stick to each")
    print("  other and polymerise, deforming the red cell into a sickle.")
    print("  One base out of 3 billion.")
    print()
    print(f"  Two names for this variant, both correct:")
    print(f"    {r['hgvs_p']:<12} HGVS -- counts from the initiator Met, which")
    print(f"                 is residue 1. This is what ClinVar records.")
    print(f"    {r['legacy_p']:<12} mature numbering -- counts after that Met is")
    print(f"                 cleaved off the finished protein. The historical")
    print(f"                 literature name, and why everyone says 'E6V'.")
    print("  HGVS = Human Genome Variation Society, who maintain the standard.")

    print()
    print("  For contrast, a third-position change in the same codon:")
    r2 = classify_snv(seq, 20, "G", "A")                   # GAG -> GAA
    print(f"  change     {r2['hgvs_c']}   (GAG -> GAA)   effect: {r2['kind']}"
          f"  {r2['hgvs_p']}")
    print("  Same codon, different position, zero consequence. That asymmetry")
    print("  is the degeneracy of the genetic code, and it is why variant")
    print("  annotation is the first filter in every clinical pipeline.")

    print()
    print(rule)
    print("5b. 'The same codon' is a claim about the ANNOTATION")
    print(rule)
    print("  Everything above assumed the coding sequence starts at index 0.")
    print("  Shift that by one or two bases and every codon boundary moves,")
    print("  so the very same base lands in a different codon at a different")
    print("  position within it. Same A>T, read three ways:")
    print()
    print(f"  {'CDS starts':<12}{'codon':>7}{'pos':>5}{'':4}{'mutant':<8}"
          f"{'effect':<12}change")
    for row in demo_frame_dependence(seq, 19, "T"):
        print(f"  index {row['offset']:<6}{row['wt_codon']:>7}{row['codon_pos']:>5}"
              f" -> {row['mut_codon']:<8}{row['kind']:<12}"
              f"{AA_NAMES[row['wt_aa']]} -> {AA_NAMES[row['mut_aa']]}")
    print()
    print("  Not a milder or stronger version of the same finding -- three")
    print("  different findings, one of which is not a missense at all but a")
    print("  lost stop codon.")
    print()
    print("  In the cell the frame is NOT a free choice: the ribosome starts at")
    print("  one particular AUG, and splicing joins the exons one particular")
    print("  way. There is a single true frame per transcript. What varies is")
    print("  our KNOWLEDGE of it -- the annotation. And that is genuinely")
    print("  plural:")
    print()
    print("   - One gene has many transcripts. Alternative splicing puts the")
    print("     same base in different codons in different isoforms, so VEP")
    print("     and SnpEff report one consequence PER TRANSCRIPT. A variant")
    print("     being missense in one and synonymous in another is routine,")
    print("     and both are true.")
    print("   - So a bare 'c.20A>T' is not interpretable. Real HGVS carries the")
    print("     transcript and its version: NM_000518.5:c.20A>T. The output")
    print("     above prints the short form only because HBB is simple.")
    print("   - MANE Select exists for exactly this reason: NCBI and EMBL-EBI")
    print("     agreeing on one canonical transcript per gene so that clinical")
    print("     labs stop talking past each other.")
    print("   - Annotation versions drift. The same variant can change")
    print("     consequence between Ensembl releases. This is a real source of")
    print("     clinical discordance, not a hypothetical one.")

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
