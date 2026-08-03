"""Population structure from a genotype matrix, with scikit-allel.

2,504 people from 26 populations, ~5Mb of chromosome 22. No labels are given
to the algorithm -- we colour by population afterwards, to check.

Run: python 05-popgen/pca.py   (after fetch_data.sh)
"""

import os
import collections
import numpy as np
import allel
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
RULE = "-" * 74

SUPERPOP_COLOR = {
    "AFR": "#c44e52", "EUR": "#4c72b0", "SAS": "#dd8452",
    "EAS": "#55a868", "AMR": "#8172b3",
}
SUPERPOP_NAME = {
    "AFR": "African", "EUR": "European", "SAS": "South Asian",
    "EAS": "East Asian", "AMR": "Admixed American",
}


def load_panel():
    """sample -> (population, superpopulation)"""
    panel = {}
    with open(os.path.join(DATA, "samples.panel")) as f:
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                panel[parts[0]] = (parts[1], parts[2])
    return panel


def ld_prune(gn, size=500, step=200, threshold=0.1, n_iter=5):
    """Thin to roughly independent markers.

    Nearby variants are correlated (linkage disequilibrium). Left in, one large
    LD block dominates PC1 and you get an axis describing a chunk of one
    chromosome instead of ancestry. This is multicollinearity, and the failure
    mode looks plausible rather than broken -- which is why it catches people.
    """
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step,
                                             threshold=threshold)
        n_before, n_after = gn.shape[0], np.count_nonzero(loc_unlinked)
        gn = gn.compress(loc_unlinked, axis=0)
        print(f"    iteration {i + 1}: {n_before:,} -> {n_after:,} variants")
        if n_after == n_before:
            break
    return gn


def main():
    print(RULE)
    print("1. Load genotypes")
    print(RULE)
    callset = allel.read_vcf(os.path.join(DATA, "chr22.slice.vcf.gz"),
                             fields=["samples", "calldata/GT", "variants/POS"])
    samples = list(callset["samples"])
    gt = allel.GenotypeArray(callset["calldata/GT"])
    print(f"  genotype array: {gt.shape[0]:,} variants x {gt.shape[1]:,} samples")
    print(f"  each cell holds two alleles (you have two copies of chr22)")
    print(f"  first sample, first 5 variants:\n    {gt[:5, 0].tolist()}")

    panel = load_panel()
    pops = [panel.get(s, ("?", "?"))[0] for s in samples]
    supers = [panel.get(s, ("?", "?"))[1] for s in samples]
    counts = collections.Counter(supers)
    print()
    for sp, n in counts.most_common():
        print(f"  {sp} {SUPERPOP_NAME.get(sp, ''):<18} {n:>5} samples")

    print()
    print(RULE)
    print("2. Genotypes -> a number matrix")
    print(RULE)
    # 0/0 -> 0, 0/1 -> 1, 1/1 -> 2. Count of alternate alleles.
    gn = gt.to_n_alt(fill=0)
    print(f"  alt-allele count matrix: {gn.shape}")
    print(f"  first sample, first 20 variants: {gn[:20, 0].tolist()}")
    print("  This is the whole trick -- from here it is only linear algebra.")

    print()
    print(RULE)
    print("3. Filter: common, non-singleton variants")
    print(RULE)
    ac = gt.count_alleles()
    # Already MAF-filtered at download, but re-derive it here so the step is
    # visible rather than hidden in a shell script.
    flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 0.05 * ac.sum(axis=1).max())
    gn = gn.compress(flt, axis=0)
    print(f"  biallelic with MAF > 5%: {gn.shape[0]:,} variants retained")
    print("  Rare variants carry no information about shared structure; they")
    print("  only add noise, and they are the ones most likely to be errors.")

    print()
    print(RULE)
    print("4. LD pruning")
    print(RULE)
    gn = ld_prune(gn)
    print(f"  final: {gn.shape[0]:,} roughly independent markers")

    print()
    print(RULE)
    print("5. PCA")
    print(RULE)
    # Patterson scaling: divide by sqrt(p(1-p)), weighting each variant by its
    # expected drift variance rather than treating all variants equally.
    coords, model = allel.pca(gn, n_components=10, scaler="patterson")
    var = model.explained_variance_ratio_ * 100
    print(f"  coords: {coords.shape}  (one row per person, one column per PC)")
    print()
    print("  variance explained:")
    for i in range(6):
        print(f"    PC{i + 1}  {var[i]:>5.2f}%  {'#' * int(var[i] * 4)}")
    print(f"    (top 10 PCs together: {var[:10].sum():.1f}%)")
    print()
    print("  Note how small these are. ~85-90% of human genetic variation is")
    print("  WITHIN populations, not between them. PCA finds the structured")
    print("  minority -- the picture is real, but it is a picture of a few")
    print("  percent of the variance.")

    print()
    print(RULE)
    print("6. Did geography fall out? (mean position per superpopulation)")
    print(RULE)
    supers_arr = np.array(supers)
    print(f"  {'group':<20}{'PC1':>9}{'PC2':>9}{'PC3':>9}")
    order = sorted(counts, key=lambda s: np.mean(coords[supers_arr == s, 0]))
    for sp in order:
        m = supers_arr == sp
        print(f"  {SUPERPOP_NAME.get(sp, sp):<20}"
              f"{coords[m, 0].mean():>9.1f}{coords[m, 1].mean():>9.1f}"
              f"{coords[m, 2].mean():>9.1f}")
    print()
    print("  PC1 should isolate AFR: African populations retain the most")
    print("  variation, because everyone else descends from a subset that")
    print("  left -- a bottleneck that discarded diversity.")
    print("  PC2 should separate EAS from EUR, with SAS lying between them.")

    print()
    print(RULE)
    print("7. Within South Asia")
    print(RULE)
    pops_arr = np.array(pops)
    sas_names = {"GIH": "Gujarati (Houston)", "PJL": "Punjabi (Lahore)",
                 "BEB": "Bengali (Bangladesh)", "STU": "Sri Lankan Tamil",
                 "ITU": "Indian Telugu (UK)"}
    sas = [p for p in sas_names if (pops_arr == p).any()]
    if sas:
        print(f"  {'population':<24}{'n':>4}{'PC1':>9}{'PC2':>9}{'PC3':>9}{'PC4':>9}")
        for p in sorted(sas, key=lambda p: coords[pops_arr == p, 2].mean()):
            m = pops_arr == p
            print(f"  {sas_names[p]:<24}{m.sum():>4}"
                  f"{coords[m, 0].mean():>9.1f}{coords[m, 1].mean():>9.1f}"
                  f"{coords[m, 2].mean():>9.1f}{coords[m, 3].mean():>9.1f}")
        print()
        print("  These five are a cline, not five clusters -- South Asian")
        print("  ancestry is a gradient of mixture between an ancestral North")
        print("  Indian component (related to West Eurasians) and an ancestral")
        print("  South Indian one. Reich's book is the readable account.")

    print()
    print(RULE)
    print("8. Plot")
    print(RULE)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, (x, y) in zip(axes, [(0, 1), (2, 3)]):
        for sp in ["AFR", "EUR", "SAS", "EAS", "AMR"]:
            m = supers_arr == sp
            if not m.any():
                continue
            ax.scatter(coords[m, x], coords[m, y], s=9, alpha=0.65,
                       c=SUPERPOP_COLOR[sp], label=SUPERPOP_NAME[sp],
                       edgecolors="none")
        ax.set_xlabel(f"PC{x + 1}  ({var[x]:.2f}%)")
        ax.set_ylabel(f"PC{y + 1}  ({var[y]:.2f}%)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("Global structure", loc="left")
    axes[1].set_title("Finer structure", loc="left")
    axes[0].legend(frameon=False, markerscale=2, fontsize=9)
    fig.suptitle("1000 Genomes, chr22 — population structure from an unlabelled "
                 "genotype matrix", x=0.02, ha="left", fontsize=12)
    fig.tight_layout()
    out = os.path.join(DATA, "pca.png")
    fig.savefig(out, dpi=140)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
