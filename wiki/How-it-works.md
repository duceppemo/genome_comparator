# How it works

## The pipeline
1. **Sketch**: each genome (assembly, or all the reads of a sample) is reduced by `mash sketch` to a *sketch*: the
   smallest hash values of all its k-mers (sequences of length k, 21 by default). 10,000 hashes are kept per genome by
   default (`--sketch-size`). Sketches are small and fast to compare.
2. **Distances**: `mash triangle` compares every pair of sketches. The fraction of shared hashes estimates the
   fraction of shared k-mers (the Jaccard index), which Mash converts into a *Mash distance*: an estimate of the
   mutation rate between the two genomes.
3. **Trees**: trees are built from the distance matrix (UPGMA always; NJ and ME optionally), and optionally
   supported by bootstrap replicates.
4. **PCoA**: the distance matrix is projected on a few axes for an overview of the structure of the dataset.

## Interpreting Mash distances
The Mash distance is roughly **1 − ANI** (average nucleotide identity) for related genomes:

| Mash distance | ≈ ANI | Typical meaning (bacteria) |
|---|---|---|
| 0 – 0.001 | > 99.9% | Same strain or very closely related isolates |
| 0.001 – 0.05 | 95 – 99.9% | Same species (e.g. the two main lineages of *L. monocytogenes* are ~0.04 apart) |
| ~0.05 | ~95% | Usual species boundary |
| > 0.1 | < 90% | Different species; distances become less precise |

* Mash distances are estimates. They are most reliable for closely to moderately related genomes; large distances
  are less precise and should not be over-interpreted. For accurate ANI values, use a dedicated tool such as
  [skani](https://github.com/bluenote-1577/skani) or [FastANI](https://github.com/ParBLiSS/FastANI)
  (see [Related tools](Related-tools)).
* Genome completeness matters: an incomplete or contaminated assembly, or reads with low coverage or many
  sequencing errors, look more distant than they are. Check `sample_stats.tsv` for unusual sizes or numbers of
  contigs.

### The p-value warning
For every pair of genomes, Mash computes the probability of seeing that many shared hashes by chance, given the
genome sizes. The log reports the largest p-value of the dataset and warns above 0.01: some pairs share no more
k-mers than unrelated genomes would. This usually means that a sample is a different organism (wrong species,
contamination, mislabelled file), or that the sketches are too small for very large genomes.

## Resolution: what Mash cannot see
A sketch holds a sample of the k-mers of a genome. Two isolates that differ by a handful of SNPs have almost exactly
the same k-mers, and their sketches may be identical or differ by only a few hashes. Their distance is then 0 or
dominated by chance, and the branching order between them in a tree is not reliable (bootstrap support will be low).
A larger `--sketch-size` improves the resolution, at the cost of time and disk space.

Mash is ideal to compare hundreds or thousands of genomes, check species identity, spot outliers or contamination,
and see the overall structure of a dataset. For outbreak investigations, where isolates differ by a few SNPs, use a
SNP-based pipeline (see [Related tools](Related-tools)).

## Choosing a tree
| Tree | File | Rooted | Assumptions | Speed |
|---|---|---|---|---|
| UPGMA (hierarchical clustering) | `_hc.nwk` | Yes | Constant rate of evolution (molecular clock) | Fastest |
| Neighbour joining | `_nj.nwk` | No | None on rates | Fast; slower than ME on very large datasets |
| Minimum evolution | `_me.nwk` | No | None on rates | Fast |

* UPGMA is a clustering method: it groups genomes by similarity and is great to find clusters, but branch lengths
  and deep branching are only reliable if genomes evolve at similar rates.
* NJ and ME trees do not assume a molecular clock and are usually preferred for phylogenetic interpretation. They are
  unrooted: root them on an outgroup (a genome from a related species) or on the midpoint in your tree viewer.
* Other linkage methods can be used for the `_hc` tree with `--linkage` (`ward`, `complete`, `single`, `weighted`).
  `ward` was the default before version 0.3.0; it assumes Euclidean distances, which Mash distances are not.

## Reading the PCoA
* Each axis shows the percentage of the total variation it explains. The first two axes are plotted; the first three
  are saved in `all_dist_PCoA.tsv`.
* Groups that overlap on the first two axes can be separated on the third one.
* A very distant genome (another genus, a contaminated assembly) can take a whole axis for itself and squeeze all the
  other genomes together. Check it, and remove it or analyse it separately if needed.
* Mash distances are not Euclidean, so some variation cannot be represented exactly (negative eigenvalues); this is
  expected and does not affect the plot.
