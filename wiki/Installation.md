# Installation

## Requirements
* Linux or macOS
* [Mash](https://github.com/marbl/Mash) 2.0 or later (installed automatically with conda)
* Python 3.10 or later with numpy, pandas, scipy, scikit-bio (≥ 0.7), plotly and openpyxl

## With pip
```
pip install genome-comparator
genome-comparator --version
```
Mash is not available from PyPI and must be installed separately, for example with
`conda install -c bioconda mash`, or from the [Mash releases](https://github.com/marbl/Mash/releases).
`genome-comparator` finds `mash` in your `PATH`, or next to the Python interpreter it runs with.

## From source with conda
```
git clone https://github.com/duceppemo/genome_comparator
cd genome_comparator
conda env create -f environment.yml
conda activate genome_comparator
pip install --no-deps .
genome-comparator --version
```

`environment.yml` only uses the `conda-forge` and `bioconda` channels. It requires `gsl>=2.8` because older
bioconda builds of Mash fail to start with newer GSL versions (see [Troubleshooting](Troubleshooting)).

## Development install
Use an editable install so changes to the code are used right away, then run the tests:
```
pip install --no-deps -e .
pytest
```
The end-to-end tests are skipped if `mash` is not in your `PATH`.

## Without installing
The main command can also be run directly from the cloned folder, as long as the dependencies are available:
```
python3 -m genome_comparator -h
```
The other tools need an installed package. From a clone, an editable install (`pip install --no-deps -e .`) is the
closest to running from the folder: `git pull` updates the commands without reinstalling.

## Updating from 0.4.7 or earlier
The scripts at the root of the repository were removed: run the installed commands instead, with the same options.

| Removed script | Command |
|---|---|
| `python3 mash_phylo.py` | `genome-comparator` |
| `python3 dendrogram_from_distance_matrix.py` | `dendrogram-from-matrix` |
| `python3 tree_collapser.py` | `tree-collapser` |
| `python3 tree_renamer.py` | `tree-renamer` |

## Updating from 0.4.3 or earlier
The main command was renamed from `mash-phylo` to `genome-comparator` in version 0.4.4, with the same options.
`mash-phylo` still works but prints a deprecation warning; update your scripts.

## Updating from 0.2 or earlier
Version 0.3.0 requires Python ≥ 3.10 and scikit-bio ≥ 0.7. Recreate your environment:
```
conda env remove -n genome_comparator
conda env create -f environment.yml
```
