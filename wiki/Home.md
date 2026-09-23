<p align="center"><img src="https://raw.githubusercontent.com/duceppemo/genome_comparator/master/assets/logo.png" width="180" alt="genome_comparator logo"></p>

# genome_comparator

**genome_comparator** quickly compares and visualizes distances between organisms using whole genome information,
from assemblies (fasta) or raw sequencing reads (fastq). It uses [Mash](https://github.com/marbl/Mash) to compute
pairwise distances, then builds trees (UPGMA, neighbour joining, minimum evolution), optionally with bootstrap
support, and an interactive PCoA plot.

It scales to thousands of genomes: sketches are computed in parallel, reused between runs, and all the pairwise
distances are computed in a single multithreaded Mash call.

## Pages
* [Installation](Installation)
* [Usage](Usage): input files, options and examples
* [Output files](Output-files)
* [Bootstrap support](Bootstrap-support)
* [Other tools](Other-tools): `dendrogram-from-matrix`, `tree-collapser`, `tree-renamer`
* [Troubleshooting](Troubleshooting)
* [Changelog](Changelog)
* [Contributing](Contributing)
