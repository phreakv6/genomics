# 01 — The Central Dogma

The whole of module 1 is one sentence: **DNA → RNA → protein**. Everything below is the detail
that makes that sentence operational, i.e. lets you write code against it.

## The molecules

**DNA** is a string over the alphabet `{A, C, G, T}` — adenine, cytosine, guanine, thymine. It is
double-stranded, and the two strands are not independent: `A` pairs with `T`, `C` pairs with `G`.
So the second strand is fully determined by the first. That redundancy is not decoration; it's the
copying mechanism and the error-correction mechanism.

The strands run in opposite directions. Each has a chemical polarity, written 5′→3′ ("five prime to
three prime"). By convention you always write a sequence 5′→3′. So the partner of a strand is not
just the complement — it's the complement *read backwards*. That's the **reverse complement**, and
it is the single most-used operation in all of bioinformatics:

```
5'-ACGTTGCA-3'   (given)
3'-TGCAACGT-5'   (complement, written in its own reading direction it becomes:)
5'-TGCAACGT-3'   (reverse complement)
```

Why you care as a programmer: a read coming off a sequencer might match the reference forwards or
backwards. Half your alignments are to the reverse strand. Get this wrong and everything downstream
is silently wrong.

**RNA** is the same thing with two changes: the sugar backbone differs (irrelevant to us), and
`T` (thymine) is replaced by `U` (uracil). **Transcription** — DNA→RNA — is therefore, in code,
`s.replace("T", "U")`. That is genuinely all it is at the string level.

**Protein** is a string over a 20-letter alphabet (the amino acids). **Translation** — RNA→protein —
reads the RNA three letters at a time. Each triplet is a **codon**, and a fixed lookup table maps
each of the 4³ = 64 codons to one of 20 amino acids or to "stop".

64 codons → 21 outcomes means the code is **degenerate**: multiple codons encode the same amino
acid, usually differing in the third position. This has a large consequence you'll meet constantly
in variant analysis: a DNA change in the third position of a codon often changes nothing about the
protein. Those are **synonymous** variants. Changes that do alter the amino acid are
**non-synonymous** (or **missense**); ones that create a premature stop are **nonsense**. Most of
variant interpretation begins with sorting mutations into these buckets.

## Genes, and why real genomes are harder than the string version

A **gene** is a stretch of DNA that gets transcribed. In humans, genes are interrupted: coding
chunks (**exons**) alternate with non-coding chunks (**introns**). The full transcript is made, then
the introns are cut out and the exons stitched together — **splicing** — before translation. And the
same gene can be spliced multiple ways (**alternative splicing**), producing different proteins from
one locus.

So the naive picture "read the genome, translate it" fails on real human DNA. You need annotation —
a separate file saying where the exons are (that's what GFF/GTF files are, module 3). This is why
genomics is a data-integration problem and not just a string problem, and it's the first place your
plumbing instincts start to pay.

Also worth internalising early: only ~1–2% of the human genome codes for protein. The rest is
regulatory sequence, repeats, and things we understand poorly. "Junk DNA" was a bad name for
"we hadn't read it yet."

## Reading frames

Because translation reads in triplets, *where you start* changes everything. Starting at position
0, 1, or 2 gives three completely different proteins. Those are **reading frames**. And since the
molecule is double-stranded, the reverse complement has its own three frames. **Six reading frames**
total for any stretch of DNA. When you don't know where a gene starts, you translate all six and
look for long stretches without a stop codon — an **open reading frame** (ORF).

Translation starts at `AUG` (which codes methionine, and doubles as the start signal) and runs until
it hits `UAA`, `UAG`, or `UGA` — the three stop codons.

Two practical consequences that bite immediately:

**Frames have different lengths.** Starting 1 or 2 bases in leaves a trailing 1–2 bases that can't
form a codon. They get discarded — two bases are not a codon and nothing can be said about them. A
93bp sequence yields 31, 30, and 30 residues in frames +1, +2, +3, not 31 three times.

**The reverse frames are not positionally paired with the forward ones.** Frame −1 starts at offset
0 of the *reverse complement*, which is the far end of the input. So index 5 in frame +1 and index 5
in frame −1 sit at opposite ends of the molecule. Naming them ±1/2/3 makes them look like aligned
columns; they are six independent readings. When you find something on the minus strand, its offset
in the reverse-complemented string you searched must be mapped back with `L - end` to be a position
on the original sequence. Skip that and the feature silently lands at the wrong end — the
reading-frame version of the coordinate bugs in module 03. `find_orfs()` handles this and proves it
by round-tripping: slice the original at the reported coordinates, reverse-complement, translate,
and check you get the peptide back.

## Variation — the vocabulary you need for modules 3–5

- **SNP / SNV** — single nucleotide polymorphism/variant. One base differs. `A` where the reference
  says `G`. The overwhelming majority of variants.
- **Indel** — a small insertion or deletion. If its length isn't a multiple of 3 inside a coding
  region it causes a **frameshift**, which garbles every downstream codon — usually catastrophic.
- **Structural variant (SV)** — large rearrangements: deletions, duplications, inversions,
  translocations, typically >50bp. Hardest to detect, and where long reads (Nanopore/PacBio) beat
  short reads (Illumina).
- **Reference genome** — an agreed-upon consensus sequence that everyone reports variants relative
  to. It is not "a normal person's genome"; it's a mosaic, and it has errors and biases.
  A "variant" is always *a difference from the reference*, never an absolute fact.
- **Allele** — one of the alternative versions at a position. You have two copies of each
  chromosome, so at any position you're **homozygous** (both copies same) or **heterozygous**
  (different). That 0/1 vs 1/1 distinction is literally what a VCF genotype field records.

## What to read

- Mukherjee, *The Gene* — the stretch on Watson/Crick through the cracking of the genetic code.
  Read for the story of how this was figured out; it makes the facts stick.
- Alberts, *Molecular Biology of the Cell*, Ch. 4 (DNA structure), Ch. 6 ("How Cells Read the
  Genome: From DNA to Protein" — this is the central dogma proper), Ch. 7 (regulation).
  Reference, not cover-to-cover. ~80–100 pages skimmed.

## Then

`central_dogma.py` implements every operation above from scratch — complement, reverse complement,
transcription, the codon table, translation, six-frame translation, ORF finding — and runs them on a
real human gene fragment. Read the code alongside this document; each function maps to a section here.
