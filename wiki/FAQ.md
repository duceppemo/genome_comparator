# FAQ

### How do I add new genomes to an existing analysis?
Put the new files in the input folder and run the same command with the same output folder. Existing sketches are
reused, so only the new genomes are sketched. Distances and trees are recomputed for all genomes.

### Can I mix assemblies and reads?
Yes, as long as each sample is either assemblies or reads. Read sketches can include sequencing errors: use
`--min-copies` (default 2) to filter k-mers seen only once, and expect distances involving low-coverage samples to be
slightly larger.

### How are sample names chosen?
From the file names, without the extension; read suffixes such as `_R1`/`_R2` are also removed from fastq files.
See [Usage](Usage#input-files).

### Why is a sample missing from the results?
Check `sample_stats.tsv`: a `failed` sample could not be read by Mash (empty, truncated or not a fasta/fastq file).
Also check the log for "Ambiguous sample names".

### How do I root the NJ or ME tree?
They are unrooted. Root them in your tree viewer on an outgroup (a genome from a related species) or on the midpoint.
The UPGMA tree (`_hc.nwk`) is already rooted.

### Why is the bootstrap support low?
Usually because the genomes are very closely related: their sketches differ by very few hashes. See
[Bootstrap support](Bootstrap-support#interpreting-support-values).

### Can I get ANI values?
ANI is roughly 1 − Mash distance for related genomes. For accurate ANI values, use skani or FastANI
(see [Related tools](Related-tools)).

### Which k-mer size and sketch size should I use?
The defaults (k = 21, 10,000 hashes) work well for bacterial genomes. A larger sketch gives more precise distances
between very close genomes. See [Usage](Usage#choosing-parameters).

### Can I use a distance matrix from another tool?
Yes: `dendrogram-from-matrix` builds the same trees and PCoA from any square distance matrix
(see [Other tools](Other-tools)).

### Can I colour the PCoA by my own groups?
Yes: `--pcoa --metadata metadata.tsv --color-by COLUMN`. The first column of the metadata file must hold the sample
names. Each group gets its own colourblind-friendly colour and marker shape.

### Can I use it on something other than bacteria?
Mash works with any genomes, but genome_comparator has mostly been used and tested on bacterial genomes.
