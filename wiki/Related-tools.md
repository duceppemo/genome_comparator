# Related tools

genome_comparator is a quick way to get distances, trees and an overview of a set of genomes. Depending on the
question, these tools may be a better fit or a good complement:

| Question | Tool | Notes |
|---|---|---|
| Trees from Mash distances | [mashtree](https://github.com/lskatz/mashtree) | Similar approach; NJ trees, with a bootstrap script |
| Accurate ANI values | [skani](https://github.com/bluenote-1577/skani), [FastANI](https://github.com/ParBLiSS/FastANI) | More precise than Mash distances, especially for distant genomes |
| Remove redundant genomes | [dRep](https://github.com/MrOlm/drep) | Dereplication with Mash, then ANI |
| Search large databases, metagenomes | [sourmash](https://github.com/sourmash-bio/sourmash) | FracMinHash sketches; search, compare and taxonomic classification |
| Strain clusters for surveillance | [PopPUNK](https://github.com/bacpop/PopPUNK) | Clusters bacterial genomes using core and accessory distances |
| Taxonomic classification | [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | Places genomes in the GTDB taxonomy |
| Outbreak resolution (few SNPs) | [Snippy](https://github.com/tseemann/snippy), [vSNP](https://github.com/USDA-VS/vSNP) | SNP-level comparison against a reference |
| Distance matrix to tree, very large datasets | [rapidNJ](https://github.com/somme89/rapidNJ), [FastME](http://www.atgc-montpellier.fr/fastme/) | Use with `mash-phylo --phylip` |

## When to use genome_comparator
* A fast first look at hundreds or thousands of genomes: structure, clusters and outliers.
* Checking species identity and spotting mislabelled or contaminated samples before a detailed analysis.
* Mixed inputs: assemblies and raw reads in the same run.
* Interactive exploration of the results with metadata (PCoA).

For the limits of Mash distances, see [How it works](How-it-works#resolution-what-mash-cannot-see).
