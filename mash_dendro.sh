#!/bin/bash


##################################
#                                #
#   Make tree from Mash output   #
#                                #
##################################


# The "export" statement makes the variables available to subprocesses,
# required to use them with GNU parallel
export baseDir=""${HOME}"/analyses/colletotrichum"  # analysis folder
export read_folder=/media/30tb_raid10/data/colletotrichum/genbank_assemblies  # make sure it's populated before running the rest
export maxProc=48  # Set maximum number of parralel files processed
export cpu=$(nproc)
[ "$maxProc" -gt "$cpu" ] && export maxProc="$cpu"  # You were too greedy!

# Create output folder for mash sketches
[ -d "${baseDir}"/mash/sketches ] || mkdir -p "${baseDir}"/mash/sketches

# Create a sketch for each file
# Creating sketches is not required, but make the distance calculations much faster
function sketch_it()
{
    sample=$(cut -d "." -f 1 <<< $(basename "$1"))

    mash sketch \
        -k 21 \
        -s 10000 \
        -p $((cpu/maxProc)) \
        -o "${baseDir}"/mash/sketches/"${sample}".msh \
        "$1"
}

export -f sketch_it

# parallel process all the files
find "$read_folder" -type f -name "*.fna" |
    parallel    --bar \
                --env sketch_it \
                --env baseDir \
                --env cpu \
                'sketch_it {}'

# create a list sketch files, one file per line
find "${baseDir}"/mash/sketches -type f -name "*.msh" > "${baseDir}"/mash/sketch.list

# paste all the sketches together
mash paste \
    "${baseDir}"/mash/all.msh \
    -l "${baseDir}"/mash/sketch.list

# Create output folder for the individual distances
[ -d "${baseDir}"/mash/distances ] || mkdir -p "${baseDir}"/mash/distances

# Create distance for each file
function dist_it()
{
    sample=$(cut -d "." -f 1 <<< $(basename "$1"))

    mash dist \
        -p $((cpu/maxProc)) \
        -t \
        "${baseDir}"/mash/all.msh \
        "$1" \
        > "${baseDir}"/mash/distances/"${sample}"_dist.tsv

    # remove paths in file names
    cat "${baseDir}"/mash/distances/"${sample}"_dist.tsv \
        | sed -e "s%$read_folder/%%g" -e 's/.fna//g' -e 's/#query//'\
        > "${baseDir}"/mash/distances/"${sample}"_dist.tsv.tmp

    mv "${baseDir}"/mash/distances/"${sample}"_dist.tsv.tmp \
        "${baseDir}"/mash/distances/"${sample}"_dist.tsv
}

export -f dist_it

find "$read_folder" -type f -name "*.fna" |
    parallel    --bar \
                --env dist_it \
                --env baseDir \
                --env read_folder \
                --env cpu \
                'dist_it {}'

# Create formated distance matrix for all files
header=$(head -n 1 $(find "${baseDir}"/mash/distances -type f -name "*_dist.tsv" | head -n 1))

# Concaternate distance files
echo "$header" > "${baseDir}"/mash/all_dist.tsv

for i in $(find "${baseDir}"/mash/distances -type f -name "*_dist.tsv"); do
    cat "$i" \
        | sed -e '1d' \
        >> "${baseDir}"/mash/all_dist.tsv
done

# Create NJ tree from matrix and save image and tree file
Rscript --vanilla \
    /home/bioinfo/scripts/NJ_tree_from_distance_matrix.R \
    "${baseDir}"/mash/all_dist.tsv
