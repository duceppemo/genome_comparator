<p align="center">
  <img src="assets/logo.png" width="180" alt="genome_comparator logo">
</p>

<h1 align="center">genome_comparator</h1>

<p align="center">
  <a href="https://github.com/duceppemo/genome_comparator/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/genome_comparator/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://github.com/duceppemo/genome_comparator/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/genome_comparator" alt="Release"></a>
  <a href="https://codecov.io/gh/duceppemo/genome_comparator"><img src="https://codecov.io/gh/duceppemo/genome_comparator/graph/badge.svg" alt="Coverage"></a>
  <img src="https://img.shields.io/badge/python-3.10%E2%80%933.14-blue" alt="Python 3.10–3.14">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/genome_comparator" alt="License"></a>
  <a href="https://github.com/duceppemo/genome_comparator/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
</p>

Quickly compare and visualize distances between genomes, from assemblies (fasta) or reads (fastq), using
[Mash](https://github.com/marbl/Mash). Produces a distance matrix, trees (UPGMA, neighbour joining, minimum
evolution) with optional bootstrap support, and an interactive PCoA plot.

## Installation
```
git clone https://github.com/duceppemo/genome_comparator
cd genome_comparator
conda env create -f environment.yml
conda activate genome_comparator
pip install --no-deps .
```

## Usage
```
mash-phylo -i /input/folder/ -o /output/folder/ --nj --pcoa
```

See the **[wiki](https://github.com/duceppemo/genome_comparator/wiki)** for all options, output files, bootstrap
support, the other tools and troubleshooting.

## Author
Marc-Olivier Duceppe: marc-olivier.duceppe@inspection.gc.ca

## License
[MIT](LICENSE)
