# Chapter 17 — Geography in the Genome

The matrix is built; now the eigendecomposition. This chapter runs PCA on
the 747-marker matrix, checks whether the structure it finds — with no
labels, from 5Mb of one chromosome — matches forty thousand years of human
migration history, and then spends serious time on how *not* to over-read
the resulting picture. It closes with where this machinery goes next
(GWAS), which lands on familiar ground: it is a factor model.

## Why PCA recovers geography

PCA finds the directions of maximum variance in the cloud of 2,504 points
(one per person, in 747-dimensional marker space). The reason those
directions turn out to mean something historical:

Most genetic differences between people are inherited, and inheritance is
geographically structured — for most of history, people had children with
people born near them. Populations separated by distance or barriers
accumulated independent drift in their allele frequencies, for tens of
thousands of years. So the *covariance structure* of genotypes encodes the
history of separations and mixtures — and PCA is precisely a device for
extracting dominant covariance structure. The landmark demonstration is
Novembre et al. (2008): PCA of European genotypes reproduces the *map of
Europe* — two principal components recover a person's country of origin to
within a few hundred kilometres.

The result, on our data — labels used only to annotate afterwards:

```text
5. PCA
--------------------------------------------------------------------------
  variance explained:
    PC1   6.08%  ########################
    PC2   3.17%  ############
    PC3   1.30%  #####
    PC4   1.02%  ####
    (top 10 PCs together: 14.4%)

6. Did geography fall out? (mean position per superpopulation)
--------------------------------------------------------------------------
  group                     PC1      PC2      PC3
  African                 -15.3     -1.6      0.6
  South Asian               3.9      1.4     -6.3
  European                  4.6      9.4      1.4
  Admixed American          4.8      2.1      6.2
  East Asian                8.4    -10.0     -0.2
```

<figure>
  <img src="img/pca.png" alt="PCA scatter plots of 1000 Genomes chr22: PC1 vs PC2 showing continental structure, PC3 vs PC4 showing finer structure, coloured by superpopulation">
  <figcaption>Figure 17-1. The repo's actual output: 2,504 people, 747 markers from 5Mb of chromosome 22, no labels given to the algorithm. Colouring by superpopulation afterwards reveals that the clusters it found are continental ancestry.</figcaption>
</figure>

Reading the axes against expectation:

**PC1 isolates Africa** (−15.3, versus +3.9 to +8.4 for everyone else) —
and the reason is one of the deepest facts in human genetics. African
populations retain the most genetic variation; all non-African peoples
descend from a relatively small group that left Africa ~60,000 years ago —
a **bottleneck** that discarded diversity. PC1, the single largest axis of
human genetic variance, is essentially *distance from the out-of-Africa
event*.

**PC2 separates East Asian (−10.0) from European (+9.4), with South Asian
between them at +1.4.** That intermediate position is not blur — it's
real history. South Asian ancestry is substantially a mixture of two
ancient components: **Ancestral North Indian** (related to West
Eurasians) and **Ancestral South Indian** (not closely related to any
population outside the subcontinent). A gradient of mixture between them
places SAS, truthfully, between EUR and EAS on this axis.

**Within South Asia**, PC3 resolves the mixture gradient itself:

```text
7. Within South Asia
--------------------------------------------------------------------------
  population                 n      PC1      PC2      PC3      PC4
  Sri Lankan Tamil         102      3.7      0.9     -6.9     -4.3
  Bengali (Bangladesh)      86      4.2     -0.0     -6.8     -3.9
  Indian Telugu (UK)       102      3.8      1.2     -6.5     -4.3
  Gujarati (Houston)       103      3.8      2.0     -6.0     -4.2
  Punjabi (Lahore)          96      4.0      2.5     -5.3     -3.6
```

Sorted by PC3, the five populations order Sri Lankan Tamil → Bengali →
Telugu → Gujarati → Punjabi — monotonically, running roughly **south to
north**, which is the direction the ANI/ASI mixture proportion runs. These
are not five clusters but a **cline** — a continuous gradient — recovered
from nothing but 5Mb of chromosome 22 genotypes. (Reich's *Who We Are and
How We Got Here* is the readable account of how this history was worked
out, and is particularly strong on the South Asian story.)

That the algorithm was never told any of this is the part that should
make you grin: population history, recovered from an eigendecomposition.

## The caveats, which are not optional

The plot is seductive, and the field has learned — sometimes painfully —
exactly how it gets over-read. Four correctives to internalise:

**The clusters are shallow.** Look back at the variance-explained bars:
PC1 carries 6% of the variance; the top ten PCs together carry **14.4%**.
Roughly 85–90% of human genetic variation is *within* populations, not
between them. PCA displays the between-group component precisely because
that's the structured part — but the picture is a picture of the minority
of the variance. Two random people from one population differ, on
average, almost as much as two random people from different continents.
The plot licenses conclusions about population *history*; it licenses
almost nothing about *individuals*.

**PCs are not ancestry components.** They are orthogonal directions of
variance. A point midway between two clusters might be genuinely admixed
— or might belong to an unsampled population that merely projects there.
PCA cannot distinguish those cases; model-based methods (ADMIXTURE) and
formal admixture tests (f-statistics) exist because of exactly this limit.
Note AMR's position in the table — Admixed American samples sit mid-plot,
and *only* external knowledge tells you that's admixture rather than a
distinct ancestral population.

**The axes are properties of the sample, not of humanity.** Add 500 more
European samples and the components rotate — variance-maximising
directions chase sample composition. There is no absolute coordinate
system here; every PCA plot describes its dataset.

**Sign and scale mean nothing.** PC1 flipping between runs is not a
finding. Only relative geometry — who is near whom, what orders along an
axis — carries information.

## Where this leads: GWAS, or the factor model

One step remains to connect this machinery to medical genetics. Take the
genotype matrix, add a **phenotype** column — a trait: height, LDL
cholesterol, a disease diagnosis — and regress the trait on each variant,
one at a time. That is a **genome-wide association study (GWAS)**. Two
of its structural features are now within reach of this book's tools:

- **The multiple-testing correction.** A million regressions means a
  million chances at p < 0.05; Bonferroni at ~a million effective tests
  gives the field's famous genome-wide significance threshold,
  **p < 5×10⁻⁸**.
- **The confounder is exactly what Chapter 17 just computed.** If a trait
  differs between populations for *any* reason (diet, environment,
  chance), then every population-differentiated variant — thousands of
  them — associates with the trait spuriously. The standard fix: include
  the top principal components as covariates in every regression,
  absorbing the ancestry axes so that only within-population signal
  remains.

If you run factor models on equities, this is not an analogy — it is the
same mathematical object. The PCs are the market's style factors;
population stratification is a common-factor exposure masquerading as an
idiosyncratic signal; putting PCs in the regression is factor
neutralisation. You have been doing population genetics all along, on a
different asset class.

## The end of the path

The five parts, walked end to end: DNA as a string; the algorithm that
compares strings at scale; the formats that carry the data; a real
pipeline scored honestly against truth; and structure across thousands of
genomes. Every result along the way was computed by code you can read in
an afternoon, on public data, on a laptop — which was the founding claim
of the repo this book documents: the distance between "published field"
and "thing you can build and verify yourself" is smaller in genomics than
almost anywhere else. The gaps that remain — annotation (Chapter 12), the
heavyweight callers, long reads, GWAS itself — are marked where they
stand.

## Run it

```bash
python 05-popgen/pca.py         # sections 5-8 are this chapter; writes data/pca.png
```

## Further reading

- Novembre et al. (2008), "Genes mirror geography within Europe" — three
  pages, one famous figure.
- Reich, *Who We Are and How We Got Here* — the narrative account of what
  this chapter computed, including ANI/ASI.
- Patterson, Price & Reich (2006) — eigenanalysis formalised, including
  significance tests for "is this PC real structure or noise?"
