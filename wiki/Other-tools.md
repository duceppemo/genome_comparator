# Other tools

## `dendrogram-from-matrix`
Build the same trees and PCoA plot as `mash-phylo` from any square distance matrix.
```
dendrogram-from-matrix -i matrix.tsv -o out/ [--linkage average] [--nj] [--me] [--pcoa]
                       [--metadata metadata.tsv] [--color-by COLUMN]
```
* Input: `.tsv`, `.csv`, `.xlsx` or `.xls`. The first row and first column hold the sample names.
* The matrix must be square, symmetric, complete and non-negative. Rows and columns may be in any order.
* Output files are named after the input file, e.g. `matrix_hc.nwk`, `matrix_PCoA.html`.
* Bootstrap support is not available here, because it needs the genomes.

## `tree-collapser`
Collapse clades whose average distance to their tips is smaller than a threshold. Useful to simplify large trees with
many near-identical genomes.
```
tree-collapser -i tree.nwk -o collapsed.nwk -d 0.001
```
A collapsed clade is replaced by a single tip named `<first tip> {<other tips>}`.

## `tree-renamer`
Rename the tips of a tree using a two-column tab-separated table: current name, new name.
```
tree-renamer -i tree.nwk -o renamed.nwk -r rename.tsv
```
* Only exact matches are renamed: `S1` never changes `S10`.
* Names from the table that are not found in the tree are listed in a warning.
