# Installation

## Requirements
* Linux or macOS
* [Mash](https://github.com/marbl/Mash) 2.0 or later (installed automatically with conda)
* Python 3.10 or later with numpy, pandas, scipy, scikit-bio (≥ 0.7), plotly and openpyxl

## With conda (recommended)
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
The commands can also be run directly from the cloned folder, as long as the dependencies are available
(the main command with `python3 -m genome_comparator`):
```
python3 -m genome_comparator -h
python3 dendrogram_from_distance_matrix.py -h
python3 tree_collapser.py -h
python3 tree_renamer.py -h
```

## Updating from 0.4.3 or earlier
The main command was renamed from `mash-phylo` to `genome-comparator` in version 0.4.4, with the same options.
`mash-phylo` still works but prints a deprecation warning; update your scripts.

## Updating from 0.2 or earlier
Version 0.3.0 requires Python ≥ 3.10 and scikit-bio ≥ 0.7. Recreate your environment:
```
conda env remove -n genome_comparator
conda env create -f environment.yml
```
