# Further Reading

Everything recommended across the chapters, consolidated, with the
reading order the repo's author would suggest for someone on the same
path.

## The backbone books

- **Buffalo, *Bioinformatics Data Skills*** (O'Reilly) — the book for a
  software engineer entering the field; Parts III and IV of this book are
  essentially a guided tour of its territory. Ch. 8–11 map directly onto
  chapters 8–16 here. Start it alongside Part III.
- **Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning
  Approach*** — the algorithms (Part II's material and far beyond),
  taught by making you implement them; pairs directly with Rosalind.
- **Alberts et al., *Molecular Biology of the Cell*** — the biology
  reference. Ch. 4 (DNA structure), Ch. 6 (the central dogma proper),
  Ch. 7 (regulation). Skim as a reference, not cover-to-cover.

## The narrative shelf

- **Mukherjee, *The Gene: An Intimate History*** — the story of how all of
  Part I was figured out; makes the facts stick. The Watson/Crick through
  genetic-code stretch backs Chapters 1–2; the sickle-cell story backs
  Chapter 4.
- **Reich, *Who We Are and How We Got Here*** — ancient DNA and human
  population history; unusually good on the South Asian (ANI/ASI) story
  that Chapter 17 recovers a slice of from raw genotypes.

## The serious second pass

- **Durbin, Eddy, Krogh & Mitchison, *Biological Sequence Analysis*** —
  alignment done probabilistically (pair HMMs), the mathematically mature
  version of Chapter 6. A second-pass book, worth it.

## Papers, in chapter order

- Li & Durbin (2009), *Fast and accurate short read alignment with
  Burrows–Wheeler transform* — the BWA paper. *(Ch. 7)*
- Poplin et al. (2018), *A universal SNP and small-indel variant caller
  using deep neural networks* — DeepVariant. *(Ch. 14)*
- Zook et al. — the Genome in a Bottle benchmark papers; consensus
  methodology, honestly described. *(Ch. 13, 15)*
- Krusche et al. (2019), *Best practices for benchmarking germline
  small-variant calls* — the GA4GH standard behind hap.py-style
  comparison. *(Ch. 15)*
- Patterson, Price & Reich (2006), *Population structure and
  eigenanalysis* — the normalisation and formal treatment behind
  Chapter 16's PCA.
- Novembre et al. (2008), *Genes mirror geography within Europe* — three
  pages, one famous figure. *(Ch. 17)*

## Specs, when you need the letter of the law

- [SAM/BAM and VCF specifications](https://samtools.github.io/hts-specs/) *(Ch. 9, 10)*
- [Ensembl GFF3/GTF documentation](https://www.ensembl.org/info/website/upload/gff.html) *(Ch. 12)*
- [rosalind.info](https://rosalind.info) — the problem source behind Part II.
