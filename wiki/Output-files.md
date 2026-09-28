# Output files

| File | Description |
|---|---|
| `all_dist.tsv` | Square distance matrix, tab-separated. The first cell is `#query`, as in `mash dist -t`. |
| `all_dist.phylip` | Same matrix in relaxed phylip format (`--phylip`). Spaces in names are replaced by `_`. |
| `sample_stats.tsv` | One row per sample (see below). |
| `all.msh` | All the sketches in a single file. Can be reused with `mash dist` or `mash screen`. |
| `sketches/` | One sketch per sample (`.msh`) and how it was made (`.json`), reused on the next run. |
| `tree/all_dist_hc.nwk` | Hierarchical clustering tree (UPGMA by default). |
| `tree/all_dist_nj.nwk` | Neighbour joining tree (`--nj`). |
| `tree/all_dist_me.nwk` | Balanced minimum evolution tree (`--me`). |
| `tree/all_dist_PCoA.html` | Interactive PCoA plot (`--pcoa`). Self-contained, works offline. |
| `tree/all_dist_PCoA.tsv` | PCoA coordinates of each sample on the first 3 axes (`--pcoa`). |
| `tree/all_dist_clusters.tsv` | Cluster of each sample at each distance threshold (`--clusters`). |
| `genome_comparator.log` | Log of the run. |

## `sample_stats.tsv`
| Column | Assemblies | Reads |
|---|---|---|
| `sample` | Sample name | Sample name |
| `type` | `fasta` | `fastq` |
| `files` | Number of files | Number of files (2 for paired-end) |
| `length` | Assembly size (bp) | Genome size estimated by Mash (bp) |
| `sequences` | Number of contigs | Number of reads |
| `est_coverage` | | Coverage estimated by Mash |
| `status` | `ok` or `failed` | `ok` or `failed` |

## Trees
* Tip labels are always single-quoted so sample names can contain any character. Single quotes in names are doubled
  (`'it''s'`), as the Newick standard requires.
* Branch lengths are Mash distances. In the `_hc` tree, branch lengths are half the merge heights, so the distance
  between two tips along the tree matches the clustering distance.
* NJ and minimum evolution trees are unrooted. Negative branch lengths from NJ are set to 0.
* With `--bootstrap`, support values (0–100) are written as internal node labels (see
  [Bootstrap support](Bootstrap-support)).

View the trees in [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), [iTOL](https://itol.embl.de/) or
[Dendroscope](https://github.com/husonlab/dendroscope3).

## PCoA
The principal coordinates analysis (classical multidimensional scaling) is the equivalent of a PCA for a distance
matrix. The axis titles give the percentage of variance explained. Hovering over a point shows the sample name and
the metadata columns (`--metadata`). See [How it works](How-it-works#reading-the-pcoa) to interpret it.

With `--color-by`, each category gets its own colour **and** marker shape, so groups can be told apart without
relying on colour alone:
* Colours come from the colourblind-friendly [Okabe-Ito](https://jfly.uni-koeln.de/color/) palette, in a fixed order.
  They were checked for all common colour vision deficiencies with every pair of colours side by side.
* Categories are assigned in "human" order, so a category keeps its colour and shape between runs: case is ignored
  and numbers are sorted by value (`1, 2, 10` and `st1, ST2, ST10` rather than `1, 10, 2` and `ST10, ST2, st1`).
* 6 colours × 7 shapes give 42 unique combinations. If there are more categories, the least frequent ones are grouped
  into "Other" (grey).
* Samples missing from the metadata file are shown as "Unknown" (grey open circles).

## Clusters
With `--clusters`, samples are grouped at one or more distance thresholds:
```
genome-comparator -i assemblies/ -o results/ --clusters 0.001 0.01 0.05
```
```
sample	cluster_0.001	cluster_0.01	cluster_0.05
S1	1	1	1
S2	1	1	1
S3	2	1	1
S4	3	2	1
```
* **Single linkage**: two samples share a cluster if they are linked by a chain of samples, each at a distance of at
  most the threshold from the next. Two samples of the same cluster can therefore be further apart than the threshold.
  Unlike cutting a tree, the result does not depend on the tree method or on the order of the samples.
* Clusters are numbered from 1 by decreasing size (ties: alphabetical order of their first sample), so the numbers
  are stable between runs on the same genomes. They can change when genomes are added.
* A sample with no other sample within the threshold is a cluster on its own.
* The log gives the number of clusters at each threshold and the size of the largest one.
* Thresholds are Mash distances, roughly 1 − ANI: 0.05 ≈ 95% ANI, the usual species boundary (see
  [How it works](How-it-works#interpreting-mash-distances)). Very small thresholds are limited by the resolution of
  Mash: isolates a few SNPs apart may have a distance of 0 (see
  [Resolution](How-it-works#resolution-what-mash-cannot-see)).

The table can be used as a metadata file to colour the PCoA by cluster, without sketching again:
```
dendrogram-from-matrix -i results/all_dist.tsv -o results/by_cluster --pcoa \
    --metadata results/tree/all_dist_clusters.tsv --color-by cluster_0.05
```
