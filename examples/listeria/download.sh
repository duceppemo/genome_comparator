#!/usr/bin/env bash
# Download the 22 Listeria RefSeq assemblies listed in metadata.tsv from NCBI into ./genomes/,
# one file per assembly, named after its accession (e.g. genomes/GCF_000196035.1.fna).
# Requires the NCBI datasets command line tool: conda install -c conda-forge ncbi-datasets-cli
set -euo pipefail

cd "$(dirname "$0")"
command -v datasets >/dev/null || { echo 'Install the NCBI datasets tool: conda install -c conda-forge ncbi-datasets-cli' >&2; exit 1; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

tail -n +2 metadata.tsv | cut -f1 > "$tmp/accessions.txt"
datasets download genome accession --inputfile "$tmp/accessions.txt" --include genome --filename "$tmp/listeria.zip"
unzip -q "$tmp/listeria.zip" -d "$tmp"

mkdir -p genomes
for fasta in "$tmp"/ncbi_dataset/data/GCF_*/*.fna; do
    cp "$fasta" "genomes/$(basename "$(dirname "$fasta")").fna"
done
echo "$(ls genomes/*.fna | wc -l) genomes saved in $(pwd)/genomes"
