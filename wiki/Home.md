<p align="center"><img src="https://raw.githubusercontent.com/duceppemo/genome_comparator/master/assets/logo.png" width="180" alt="genome_comparator logo"></p>

# genome_comparator

**genome_comparator** quickly compares and visualizes distances between organisms using whole genome information,
from assemblies (fasta) or raw sequencing reads (fastq). It uses [Mash](https://github.com/marbl/Mash) to compute
pairwise distances, then builds trees (UPGMA, neighbour joining, minimum evolution), optionally with bootstrap
support, and an interactive PCoA plot.

It scales to thousands of genomes: sketches are computed in parallel, reused between runs, and all the pairwise
distances are computed in a single multithreaded Mash call.

![Neighbour joining tree of 22 Listeria genomes](https://raw.githubusercontent.com/duceppemo/genome_comparator/master/assets/tree.png)

## Getting started
* [Installation](Installation)
* [Tutorial](Tutorial): 22 public *Listeria* genomes, step by step
* [Usage](Usage): input files, options and examples

## Understanding the results
* [Output files](Output-files)
* [How it works](How-it-works): interpreting distances, trees and the PCoA
* [Bootstrap support](Bootstrap-support)
* [Performance](Performance)

## Reference
* [Other tools](Other-tools): `dendrogram-from-matrix`, `tree-collapser`, `tree-renamer`
* [Related tools](Related-tools)
* [Troubleshooting](Troubleshooting)
* [FAQ](FAQ)
* [Changelog](Changelog)
* [Contributing](Contributing)

## Citing
Please cite genome_comparator: Duceppe M-O. genome_comparator: fast comparison and visualization of genome distances
with Mash. Zenodo. https://doi.org/10.5281/zenodo.22920856 (all versions; each release also has its own DOI on
Zenodo). Please also cite Mash: Ondov BD *et al.* Mash: fast genome and metagenome distance estimation using MinHash.
*Genome Biology* 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x
