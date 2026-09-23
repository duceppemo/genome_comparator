# Troubleshooting

### `"mash" was not found`
`genome-comparator` uses the `mash` found in your `PATH` first, then the one installed next to its own Python interpreter
(the `bin/` folder of the conda environment). Calling `/path/to/envs/genome_comparator/bin/genome-comparator` without
activating the environment therefore works if Mash is installed in that environment. The log shows which `mash`
was used.

If Mash cannot be found in either place, install it in the environment: `conda install -c bioconda mash`.

### `mash: error while loading shared libraries: libgsl.so.25`
Older bioconda builds of Mash are linked to GSL 2.6. Recreate the environment from `environment.yml`, which requires
`gsl>=2.8` and a matching Mash build, and do not use the `defaults` channel:
```
conda env remove -n genome_comparator
conda env create -f environment.yml
```

### `Ambiguous sample names`
Two input files give the same sample name. Rename or remove the files listed in the message.
See [Usage](Usage#input-files) for how sample names are derived from file names.

### `Could not sketch "<sample>"`
Mash could not read the file: it may be empty, truncated or not a fasta/fastq file. The sample is flagged as `failed`
in `sample_stats.tsv` and left out; the run continues if at least 3 samples remain.

### `Some distances are not significant`
The largest Mash p-value is above 0.01. Some samples may be unrelated (e.g. a contaminant or a wrong species), or the
sketch size is too small. Check `all_dist.tsv` for samples with unusually large distances, or increase `--sketch-size`.

### The tree has low bootstrap support
See [Bootstrap support](Bootstrap-support#interpreting-support-values).

### Getting more details
Run with `-v` to log every Mash command and show the full Python traceback of errors. The log is also saved in
`<output>/genome_comparator.log`.
