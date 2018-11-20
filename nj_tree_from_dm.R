#!/usr/bin/env Rscript

# load required library
if (!require("ape")){
    install.packages("ape")
} else {
    library("ape")
}

if (!require("data.table")){
    install.packages("data.table")
} else {
    library("data.table")
}

# the matrix should look like this:
    
#              sample1   sample2   sample3
#    sample1   0         0.14077   0.20428
#    sample2   0.23131   0         0.23254
#    sample3   0.20428   0.20484   0


# Get argument from command line
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Please use the matrix file as argument", call.=FALSE)
} else {
    input_file = args[1]
}

# Input file
# input_file = "/home/bioinfo/analyses/colletotrichum/mash/all_dist.tsv"
input_path = dirname(input_file)
input_name = basename(input_file)
input_name_noExt = strsplit(input_name, "[.]")[[1]][1]  # keep name without extension

# Output files
out_tree_file = paste(input_path, "/", input_name_noExt, ".tree", sep = "")
picture_type = "png"
out_picture_file = paste(input_path, "/", input_name_noExt, ".", picture_type, sep = "")
    
# import data
m <- as.matrix(fread(input_file), rownames = 1)

# Order the matrix so column are in same order as rows
m <- m[order(rownames(m)), order(colnames(m))]

# create tree
arbol <- nj(m)

# Write tree to file
write.tree(arbol, file = out_tree_file)

# Plot to file
# Create picture file handle
png(out_picture_file)

# Make tree
plot(arbol)

# Close picture file handle
dev.off()
