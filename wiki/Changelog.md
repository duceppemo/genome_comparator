# Changelog

## Unreleased
* Releases are archived on Zenodo: https://doi.org/10.5281/zenodo.22920856 (all versions). DOI badge and citation
  added to the README and `CITATION.cff`.

## 0.4.5 (2026-09-23)
* Fixed: `--me` and `--bootstrap` crashed with exactly 3 genomes.
* Fixed: `tree-renamer` and `tree-collapser` wrote bootstrap support values as quoted labels (`'95'`), which tree
  viewers do not read as support values. They are now kept as support values.
* Fixed: reading an `.xls` matrix without the `xlrd` package gave a traceback; the error now explains how to fix it.
* Fixed: the PCoA axes showed "nan%" when all genomes were identical.
* `tree-renamer` warns when several tips end up with the same name.
* Packaging: SPDX license metadata (no more setuptools deprecation warning); README images use absolute links.

## 0.4.4 (2026-09-23)
* The main command is now `genome-comparator`, the same name as the tool, with the same options.
  `mash-phylo` still works but prints a deprecation warning. `python -m genome_comparator` also works.
* Description added to `CITATION.cff` (used by Zenodo).

## 0.4.3 (2026-09-23)
* PCoA plots use a colourblind-friendly palette (Okabe-Ito) and one marker shape per category, instead of Plotly's
  default colours. Missing metadata is shown as "Unknown"; past 42 categories, the least frequent ones are grouped
  into "Other".
* Fixed: `sample_stats.tsv` left the number of sequences empty for single-sequence assemblies (e.g. complete genomes).
* New example dataset (`examples/listeria/`, 22 public genomes) and wiki pages: Tutorial, How it works,
  Performance, Related tools and FAQ.
* README with example figures, features, citation information; `CITATION.cff`, `CONTRIBUTING.md`, code of conduct,
  security policy, and issue and pull request templates.

## 0.4.2 (2026-09-23)
* `mash-phylo` now finds Mash in its own conda environment when `mash` is not in the `PATH`, so
  `/path/to/envs/genome_comparator/bin/mash-phylo` works without activating the environment.
  The log shows which `mash` was used.

## 0.4.1 (2026-09-23)
* Fixed: `--bootstrap` could hang on Python < 3.14. Worker processes were started with `fork` (the Linux default
  before Python 3.14), which can deadlock when the parent process runs threads. They now use `forkserver`
  (or `spawn` where `forkserver` is not available).
* CI: 20-minute timeout on the test jobs, and no duplicate test runs on tag pushes.

## 0.4.0 (2026-09-23)
* New `--bootstrap N` option: support values on all the trees (UPGMA, NJ, ME), from replicates sketched with
  different hash seeds. Replicate trees are built in parallel worker processes.
  See [Bootstrap support](Bootstrap-support).
* Documentation moved to this wiki, maintained in the `wiki/` folder of the repository.

## 0.3.0 (2026-09-23)
* Package layout with installable commands and a `pyproject.toml`; `conda env create -f environment.yml` works.
* Pairwise distances are computed with a single `mash triangle` call instead of one `mash dist` per sample.
* Sketches are named after the sample, so the slow file-path-to-name substitution is gone.
* Fixed: fasta files whose name contains `_R1`/`_R2` (e.g. `Iso_R10.fasta` and `Iso_R11.fasta`) were merged into
  a single sample.
* Fixed: branch lengths of the hierarchical clustering tree were rounded to 2 decimals (most Mash distances became 0).
* Fixed: re-running in the same output folder kept the old `all.msh` because `mash paste` refuses to overwrite it.
* Fixed: failures of Mash were ignored; they are now reported.
* Fixed: relative input paths starting with `../` pointed to the wrong folder.
* Fixed: `--clean` recursively deleted every `*.list` file under the output folder; it now only removes files it created.
* Fixed: `tree-renamer` used partial matches (`S1` also renamed `S10`); `tree-collapser` wrote invalid Newick when
  names contained commas.
* The default hierarchical clustering method is now `average` (UPGMA) rather than `ward`, which assumes Euclidean
  distances. Use `--linkage ward` for the previous behaviour.
* The PCA was replaced by a PCoA (the proper ordination for a distance matrix), with optional metadata colouring.
* New `--me` tree, `--phylip` output, `--metadata`/`--color-by`, `--min-copies`, `--force` and `--version` options.
* `sample_stats.txt` is now a tab-separated table: `sample_stats.tsv`.
* `tree-collapser` no longer needs ete3.
