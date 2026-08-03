# 05 — Population structure falls out of a PCA

This is the module where your existing skills stop being *transferable* and simply *are* the skill.
Everything here is linear algebra and statistics on a matrix. The biology is in the interpretation,
not the computation.

## The data: 1000 Genomes

The 1000 Genomes Project sequenced 2,504 people (Phase 3) from 26 populations grouped into five
superpopulations: AFR (African), EUR (European), **SAS (South Asian)**, EAS (East Asian), AMR
(admixed American). The SAS group includes GIH (Gujarati), PJL (Punjabi), BEB (Bengali), STU (Sri
Lankan Tamil), ITU (Indian Telugu) — so the population you belong to is directly represented, which
makes the output personally legible in a way toy data isn't.

The data is a VCF per chromosome: rows are variant sites, columns are people, cells are genotypes.

## From genotypes to a matrix

Encode each genotype as a count of alternate alleles:

```
0/0 -> 0     homozygous reference
0/1 -> 1     heterozygous
1/1 -> 2     homozygous alternate
```

That gives an integer matrix `G` of shape (variants × samples). This is the whole trick: a
population-genetics dataset is a numeric matrix, and everything downstream is standard multivariate
statistics.

Before PCA, two preprocessing steps that are not optional:

**1. Filter to common, high-quality, biallelic SNPs.** Rare variants (present in one or two people)
carry no information about shared structure and add noise. Standard: minor allele frequency > 5%.

**2. LD pruning.** Nearby variants are inherited together in blocks — **linkage disequilibrium** —
so they are highly correlated columns. Without pruning, one large LD block (or an inversion, or the
MHC region) dominates the top principal component, and you get a PC that describes one chunk of one
chromosome rather than ancestry. This is precisely the multicollinearity problem you already know
from regression; the fix is to thin to roughly independent markers. Skipping this step is the most
common way a popgen PCA goes wrong, and the failure looks plausible rather than broken.

**3. Standardise.** Center each variant and scale it. The conventional scaling divides by
`sqrt(p(1-p))` where `p` is the allele frequency — Patterson's normalisation, which weights variants
by their expected drift variance rather than treating all equally.

## What PCA is doing here, and why it works so well

PCA finds the directions of maximum variance in the sample cloud. Since most genetic variation
between people is inherited, and inheritance is geographically structured by tens of thousands of
years of migration and isolation, **the top principal components of human genotype data recover
geography.** Not approximately — the famous Novembre et al. (2008) result showed that a PCA of
European genotypes reproduces the map of Europe, to the point where you can identify someone's
country of origin from two numbers.

The standard structure you should expect to see:

- **PC1** separates African from non-African samples. This is the deepest split, because African
  populations retain the most variation — everyone else descends from a subset that left, a
  bottleneck that discarded diversity. PC1 is essentially "distance from the out-of-Africa event."
- **PC2** separates East Asian from European, with **South Asian sitting between them**, closer to
  European. That intermediate position is real history: South Asian ancestry is substantially a
  mixture of an ancestral North Indian component (related to West Eurasians) and an ancestral South
  Indian component.
- **PC3+** resolve within-group structure — and within SAS you can start to see the Gujarati /
  Punjabi / Telugu / Tamil / Bengali samples separate along a cline.

No labels are given to the algorithm. The clustering is unsupervised; you colour by known population
*afterwards*, to check. That's the part that should make you grin: population history is recoverable
from an eigendecomposition of a genotype matrix.

## Caveats worth internalising before you over-read a plot

- **PCs are not ancestry components.** They're orthogonal axes of variance. A point midway between
  two clusters may be genuinely admixed — or may be from an unsampled population that happens to
  project there. PCA cannot distinguish these; ADMIXTURE and f-statistics (Patterson, Reich) exist
  because of that limit.
- **Sample composition determines the axes.** Add 500 more Europeans and the components rotate.
  There is no absolute coordinate system; the plot describes your dataset, not humanity.
- **Sign and scale are arbitrary.** PC1 flipping between runs means nothing.
- **The clusters are shallow.** Roughly 85–90% of human genetic variation is *within* populations,
  not between them. PCA visualises the small between-group component precisely because that's the
  structured part; the picture is real but it is a picture of the minority of the variance, and it
  does not license conclusions about individuals.

## Where this leads

The same matrix, plus a phenotype column, is a **GWAS** — regress the trait on each variant, one at a
time, and correct for a million tests (Bonferroni gives the field's famous p < 5×10⁻⁸ threshold).
And the principal components you computed here go straight in as covariates, because population
structure is the classic confounder: if a trait is more common in one population, every
population-differentiated variant associates with it spuriously. Controlling for PCs is how the
field handles that. It is the exact shape of a factor model in equities.

## What to read

- Buffalo, Ch. 8–9 — the data-wrangling side.
- Novembre et al. (2008), "Genes mirror geography within Europe" — three pages, one famous figure.
- Reich, *Who We Are and How We Got Here* — narrative, and unusually good on the South Asian story
  (the ANI/ASI mixture). Very much your kind of book.

## Then

`pca.py` downloads a slice of real 1000 Genomes chromosome 22 data, does the filtering, pruning and
standardisation described above with `scikit-allel`, and plots the result coloured by population.
