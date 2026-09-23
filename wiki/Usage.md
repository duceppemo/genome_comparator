# Usage

## Quick start
```
genome-comparator -i /input/folder/ -o /output/folder/ -t 48 --nj --pcoa
```

## Input files
* The input folder is searched recursively for `.fasta`, `.fa`, `.fna`, `.fas`, `.fastq` and `.fq` files,
  optionally gzipped (`.gz`). Other files are ignored.
* The sample name is the file name without its extension. Dots elsewhere in the name are kept
  (`E.coli_K12.v2.fasta.gz` → `E.coli_K12.v2`).
* For fastq files, read suffixes are also removed: `_R1`, `_R2`, `_1`, `_2`, `_R1_001`, `_R2_001`. All the fastq files
  of a sample (e.g. paired-end reads) are combined in a single sketch. Fasta names are never trimmed.
* Ambiguous sample names stop the run with an error instead of being silently merged:
  * two fasta files with the same sample name (e.g. `S1.fasta` and `S1.fna`),
  * a sample with both fasta and fastq files,
  * the same file name in several subfolders.
* Files that Mash cannot read are reported, flagged as `failed` in `sample_stats.tsv` and left out of the analysis.
* At least 3 samples are required.

## Options
```
  -i, --input           Folder containing the fasta or fastq files (searched recursively)
  -o, --output          Folder to hold the result files
  -t, --threads         Number of threads (default: all available CPUs)
  -k, --kmer-size       k-mer size used by Mash, 1 to 32 (default: 21)
  -s, --sketch-size     Number of min-hashes per sketch (default: 10000)
  -m, --min-copies      Reads only: minimum copies of a k-mer to be kept, filters out sequencing
                        errors (default: 2)
  -b, --bootstrap N     Add bootstrap support values from N replicates to all the trees (default: 0)
  --phylip              Also save the distance matrix in phylip format
  --force               Sketch all samples again, even if up-to-date sketches exist
  --clean               Remove the individual sketch files at the end
  -v, --verbose         Show debug messages
  --version             Show the version and exit

Trees and ordination:
  --linkage             Clustering method for the "_hc" tree: average (UPGMA, default), ward,
                        complete, single or weighted
  --nj                  Also build a neighbour joining tree
  --me                  Also build a balanced minimum evolution tree (with NNI), faster than NJ on
                        very large datasets
  --pcoa, --pca         Also run a principal coordinates analysis (PCoA) and save an interactive plot
  --metadata            Tab-separated file (first column = sample name) shown when hovering over
                        PCoA points
  --color-by            Metadata column used to colour the PCoA points
```
Run `genome-comparator -h` for the full help.

## Reusing sketches
Sketches are kept in `<output>/sketches/` and reused on the next run into the same output folder if the input
files, k-mer size, sketch size and `--min-copies` did not change. Adding a few genomes to a large dataset therefore
only sketches the new ones. Use `--force` to sketch everything again, or `--clean` to delete the sketches at the end.

## Choosing parameters
* **k-mer size**: 21 works well for bacterial genomes. Larger values are more specific but more sensitive to
  sequencing errors and divergence.
* **Sketch size**: larger sketches give more precise distances between closely related genomes, at the cost of
  time and disk space. The default of 10,000 is 10 times Mash's own default.
* **Tree method**: the `_hc` tree (UPGMA) is always built and is very fast. NJ and minimum evolution trees do not
  assume a molecular clock and are usually preferred for phylogenetic interpretation.
* The log reports the largest Mash p-value. A warning is shown above 0.01: some samples may be unrelated or the
  sketches too small for reliable distances.

## Examples
Assemblies, all trees, PCoA coloured by serovar:
```
genome-comparator -i assemblies/ -o results/ --nj --me --pcoa --metadata metadata.tsv --color-by serovar
```

Paired-end reads, keeping only k-mers seen at least 3 times:
```
genome-comparator -i reads/ -o results/ -m 3
```

Trees with 100 bootstrap replicates:
```
genome-comparator -i assemblies/ -o results/ --nj --me --bootstrap 100
```

Matrix for another tree program (e.g. [rapidNJ](https://github.com/somme89/rapidNJ)):
```
genome-comparator -i assemblies/ -o results/ --phylip
rapidnj results/all_dist.phylip -i pd > results/rapidnj.nwk
```

### Metadata file format
Tab-separated; the first column holds the sample names as they appear in the results:
```
sample	serovar	source
S1	Enteritidis	poultry
S2	Typhimurium	swine
```
