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
the metadata columns (`--metadata`). Points can be coloured by a metadata column (`--color-by`).
