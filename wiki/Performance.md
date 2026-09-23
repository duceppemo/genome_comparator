# Performance

## Real datasets
| Dataset | Threads | Normal run | With `--bootstrap 100` | Peak memory |
|---|---|---|---|---|
| 22 *Listeria* genomes, ~3 Mb ([Tutorial](Tutorial)) | 8 | 3 s | 1 min 40 s | 220 MB |
| 321 *Mycobacterium bovis* assemblies, ~4.2 Mb | 64 | 7 s | 11 min 44 s | — |

Both runs used the default k-mer size (21) and sketch size (10,000), with `--nj --me --pcoa`.

* **Sketching** takes most of the time of a normal run. It is parallel (one Mash process per genome) and sketches are
  reused on the next run into the same output folder, so adding genomes to an existing dataset only sketches the new
  ones.
* **Distances** (`mash triangle`) are multithreaded and take seconds for hundreds of genomes. The number of pairs grows
  with the square of the number of genomes.
* **Bootstrap**: each replicate re-sketches all the genomes, so N replicates take about N times the sketching time.

## Tree building and PCoA
Time to build each tree from an existing distance matrix (random distances, one process, one run each):

| Genomes | UPGMA | Minimum evolution | Neighbour joining | PCoA | Peak memory |
|---|---|---|---|---|---|
| 1,000 | 0.1 s | 2.3 s | 0.1 s | 0.1 s | 0.2 GB |
| 2,500 | 0.2 s | 1.8 s | 1.9 s | 0.1 s | 0.6 GB |
| 5,000 | 0.8 s | 8.2 s | 20 s | 0.2 s | 1.7 GB |
| 9,000 | 3.1 s | 26 s | 1 min 52 s | 0.5 s | 5.3 GB |

* Neighbour joining time grows about with the cube of the number of genomes; above a few thousand genomes,
  minimum evolution (`--me`) is faster.
* A distance matrix of n genomes takes n² × 8 bytes of memory: 650 MB for 9,000 genomes. With `--bootstrap`, each
  tree-building process (up to 8) holds its own copy.
* For very large datasets, `--phylip` saves the matrix for dedicated tree programs such as
  [rapidNJ](https://github.com/somme89/rapidNJ) or [FastME](http://www.atgc-montpellier.fr/fastme/).

## Tips
* Use `-t` to match the CPUs you have (default: all available CPUs, respecting SLURM and other CPU limits).
* Keep the output folder between runs to reuse sketches; use `--clean` only when you are done.
* Increase `--sketch-size` only if you need a finer resolution between very close genomes: time and disk usage grow
  with it.
