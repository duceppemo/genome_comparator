# genome_comparator
A tool to quickly compare and visualize distances between organisms using whole genome information from assemblies
(fasta) or raw sequening data fastq).


## Author
Marc-Olivier Duceppe: marc-olivier.duceppe@canada.ca

## Usage
A typical command to compare bacteria genomes assemblies:
```
python3 mash_phylo.py -i /input/folder/ -o /output/folder/ -t 48 --nj 
```
The program's help:
```
usage: mash_phylo.py [-h] -i /input/folder -o /output/folder [-t 48]
                     [-k KMER_SIZE] [-s SKETCH_SIZE] [--nj] [--pca]

Create dendrogram from a square distance matrix

optional arguments:
  -h, --help            show this help message and exit
  -i /input/folder, --input /input/folder
                        Input folder containing the fasta or fastq files
  -o /output/folder, --output /output/folder
                        Folder to hold the result files
  -t 48, --threads 48   Number of CPU.Default is maximum CPU available(48)
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        kmer size to use for Mash.Default 21.
  -s SKETCH_SIZE, --sketch_size SKETCH_SIZE
                        Sketch size to use for MashDefault 10,000.
  --nj                  Output neighbour joining (nj) tree. Default
                        "False".Warning: can take several hours to complete on
                        big matricese.g. it takes about 8h to run on a 9,000 x
                        9,000 matrix
  --pca                 Output PCA. Default "False".Not very useful when too
                        many samples.
                        Still experimental
```

##Output

* Square distance matrix (tab-separated file).
* Hierarichal clustering dendrogram. Created very quickly.
* Neighbour-Joining tree (Optional, but recommended). Comparing 9,000 assemblies takes about 9h to complete.

Visualize the tree/dendrogram in you favorite tree viewing tool. I like to use 
[FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or 
[dendroscope](http://ab.inf.uni-tuebingen.de/software/dendroscope/welcome.html).
