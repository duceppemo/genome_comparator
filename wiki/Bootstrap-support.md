# Bootstrap support

```
mash-phylo -i assemblies/ -o results/ --nj --me --bootstrap 100
```

## How it works
Mash distances are not computed from an alignment, so there are no alignment columns to resample. Instead, each
replicate sketches all the samples again with a different hash seed (`mash sketch -S`). This selects a different
random subset of k-mers from every genome, which gives a slightly different distance matrix. The same trees are
rebuilt from each replicate matrix. mashtree uses the same approach.

The support of a clade is the percentage of replicate trees that contain it:
* `_hc` tree (rooted): the same group of tips must form a clade in the replicate tree.
* `_nj` and `_me` trees (unrooted): the same bipartition of the tips must be present in the replicate tree.

Support values are written as internal node labels in the Newick files, e.g. `(('A':0.01,'B':0.01)95:0.02,...)`.
FigTree and iTOL can display them on the branches.

## Runtime
* Each replicate re-sketches every sample, so N replicates take roughly N times as long as the sketching step of a
  normal run.
* Sketching and distances use all the threads. Replicate trees are built in parallel worker processes (up to 8) while
  the next replicate is being sketched.
* Each tree-building process holds one distance matrix in memory: about 650 MB for 9,000 samples.
* Replicate sketches are temporary: they are deleted after each replicate and not reused between runs.

## Interpreting support values
Very closely related genomes, such as isolates from an outbreak, can differ by only a few SNPs. Their sketches then
differ by very few hashes, if any, so their clades often get low support. This reflects the resolution limit of
Mash, not a problem with the run. Use a SNP-based approach to resolve relationships at that level. Deeper clades are
usually well supported.
