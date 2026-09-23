# genome_comparator
A tool to quickly compare and visualize distances between organisms using whole genome information from assemblies
(fasta) or raw sequencing data (fastq), using [Mash](https://github.com/marbl/Mash).

## Author
Marc-Olivier Duceppe: marc-olivier.duceppe@inspection.gc.ca

## Installation
```
git clone https://github.com/duceppemo/genome_comparator
cd genome_comparator
conda env create -f environment.yml
conda activate genome_comparator
pip install --no-deps .
mash-phylo -h
```
The scripts can also be run directly from the cloned folder without `pip install`
(e.g. `python3 mash_phylo.py -h`).

## Usage
A typical command to compare bacterial genome assemblies:
```
mash-phylo -i /input/folder/ -o /output/folder/ -t 48 --nj --pcoa
```

* The input folder is searched recursively for `.fasta`, `.fa`, `.fna`, `.fas`, `.fastq` and `.fq` files
  (optionally gzipped).
* The sample name is the file name without its extension. For fastq files, read suffixes
  (`_R1`, `_R2`, `_1`, `_2`, `_R1_001`, ...) are also removed, and all the files of a sample (e.g. paired-end reads)
  are combined in a single sketch.
* Ambiguous sample names (e.g. `S1.fasta` and `S1.fna`, or a sample with both fasta and fastq files) stop the run
  with an error instead of being silently merged.
* Files that Mash cannot read are reported, flagged as `failed` in `sample_stats.tsv` and left out of the analysis.
* Sketches are kept in the output folder and reused on the next run if the input files and parameters did not change,
  so adding a few genomes to a large dataset only sketches the new ones. Use `--force` to sketch everything again.

### Main options
```
  -i, --input           Folder containing the fasta or fastq files (searched recursively)
  -o, --output          Folder to hold the result files
  -t, --threads         Number of threads (default: all available CPUs)
  -k, --kmer-size       k-mer size used by Mash (default: 21)
  -s, --sketch-size     Number of min-hashes per sketch (default: 10000)
  -m, --min-copies      Reads only: minimum copies of a k-mer to be kept (default: 2)
  --linkage             Clustering method for the "_hc" tree: average (UPGMA, default), ward, complete,
                        single or weighted
  --nj                  Also build a neighbour joining tree
  --me                  Also build a balanced minimum evolution tree
  --pcoa                Also run a PCoA and save an interactive html plot (--pca is an alias)
  --metadata            Tab-separated file (first column = sample name) shown when hovering over PCoA points
  --color-by            Metadata column used to colour the PCoA points
  --phylip              Also save the distance matrix in phylip format
  --force               Sketch all samples again
  --clean               Remove the individual sketch files at the end
```

## Output
* `all_dist.tsv`: square distance matrix (tab-separated).
* `all_dist.phylip`: same matrix in phylip format (`--phylip`), e.g. for [rapidNJ](https://github.com/somme89/rapidNJ)
  or [FastME](http://www.atgc-montpellier.fr/fastme/).
* `sample_stats.tsv`: assembly size and number of contigs, or estimated genome size, number of reads and coverage
  for reads.
* `all.msh`: all the sketches in a single file. Can be reused with `mash dist` or `mash screen`.
* `tree/all_dist_hc.nwk`: hierarchical clustering tree. Created very quickly.
* `tree/all_dist_nj.nwk` and `tree/all_dist_me.nwk`: neighbour joining and minimum evolution trees (optional).
* `tree/all_dist_PCoA.html` and `tree/all_dist_PCoA.tsv`: interactive PCoA plot and coordinates (optional).
* `genome_comparator.log`: log of the run.

Visualize the trees in your favorite tree viewer, e.g. [FigTree](http://tree.bio.ed.ac.uk/software/figtree/),
[iTOL](https://itol.embl.de/) or [Dendroscope](https://github.com/husonlab/dendroscope3).
Tip labels are always quoted in the Newick files so sample names can contain any character.

## Other tools
* `dendrogram-from-matrix -i matrix.tsv -o out/ [--nj] [--me] [--pcoa]`: build trees and a PCoA plot from any
  square distance matrix (`.tsv`, `.csv`, `.xlsx` or `.xls`).
* `tree-collapser -i tree.nwk -o collapsed.nwk -d 0.001`: collapse clades whose average distance to their tips
  is smaller than a threshold. Collapsed clades are renamed `<first tip> {<other tips>}`.
* `tree-renamer -i tree.nwk -o renamed.nwk -r rename.tsv`: rename tips using a two-column tab-separated table
  (current name, new name). Only exact matches are renamed.

## Tests
```
pip install --no-deps -e .
pytest
```

## Changes in 0.3.0
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
