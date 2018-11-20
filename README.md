# genome_comparator
Compare and visualize distances between genomes quickly

## Usage

1. run the script 'mash_dendro.sh', which will produce a square distance matrix (tab-separated file). In it's current version, you need to modify the first few lines to point to right paths. There are hardcoded paths too, which will need to be changed. This will be all cleaned up in a future version. At the end of the bash script, the R script 'nj_tree_from_dm.R' is called to make a nj tree file.

2. Optional. If you want, you can use the python script 'dendro_from _dm.py' to make the tree. This script actually makes 2 trees: (1) a dendrogram based on hierarichal clustering; and (2) a NJ tree. The idea is that the dendrogram is very quick to be computed, which allows to have a quick idea of "phylogeny" of the samples.

3. Visualize the tree/dendrogram in you favorite tree viewing tool. I like to use [FigTree](http://tree.bio.ed.ac.uk/software/figtree/).
