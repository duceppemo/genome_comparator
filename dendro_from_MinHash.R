#!/usr/bin/env Rscript

# Input: signature files from sourmash (json)
# TODO -> interface with parent script (the argument part)
# TODO -> add option to choose distance method maybe?

# load required library
if (!require("rjson")){
    install.packages("rjson", dependencies = TRUE)
    library("rjson")
} else {
    library("rjson")
}

if (!require("dplyr")){
    install.packages("dplyr", dependencies = TRUE)
    library("dplyr")
} else {
    library("dplyr")
}

if (!require("pvclust")){
    install.packages("pvclust", dependencies = TRUE)
    library("pvclust")
} else {
    library("pvclust")
}

# the matrix should look like this:
    
#              sample1   sample2   sample3
#    sample1   0         0.14077   0.20428
#    sample2   0.23131   0         0.23254
#    sample3   0.20428   0.20484   0


# Get argument from command line
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1 | ! dir.exists(input_folder)) {
    stop("Please use sourmash signature folder as input", call.=FALSE)
} else {
    input_file = args[1]
}

input_folder <- "/home/bioinfo/analyses/colletotrichum/sourmash/signatures"

# Master datafram
df <- data.frame()
file_list <- list.files(input_folder, full.names = TRUE, pattern = "*.sig")
header_list <- NULL
for (input_file in file_list){
    #parse_json(input_file)
    input_name = basename(input_file)
    input_name_noExt = strsplit(input_name, "[.]")[[1]][1]  # keep name without extension
    
    json_content<- fromJSON(file = input_file)
    tmp_df <- data.frame(json_content[[1]]$signatures[[1]]$mins)
    colnames(tmp_df) <- input_name_noExt
    #df <- merge(df, tmp_df, by.x = 0, by.y = 0, all = TRUE)
    if (is.null(header_list)) {  # The fisrt iteration
        df <- tmp_df
    } else {
        df <- merge(df, tmp_df, by = 'row.names', all = TRUE)
        df <- select(df, -1)  # remove the 'row.names' column added at the first position in the data frame
    }
    header_list <- c(header_list, input_name_noExt)
}

# sort all dataframe columns independently
df <- apply(df, 2, sort, decreasing=F)

jaccard <- function(x){
    # Create an empty similarity matrix
    jaccard.sim <- as.matrix(matrix(data = NA, nrow = ncol(x), ncol = ncol(x), dimnames = list(colnames(x), colnames(x))))
    
    for (i in 1:(ncol(x)-1)){
        jaccard.sim[i,i] = 1.0
        for (j in (i+1):ncol(x)){
            #  m11 / (m01 + m10 + m11)
            intersection.length <- length(intersect(x[,i], x[,j]))
            union.length <- length(union(x[,i], x[,j]))
            
            sim = intersection.length / union.length
            #sim.other = intersection.length / (length(x[,i]) + length(x[,j]) - intersection.length)
            #print(c(sim, sim.other))  # Gives same results
            jaccard.sim[i,j] = sim
            jaccard.sim[j,i] = sim
        }
    }
    jaccard.dist <- 1 - jaccard.sim
    attr(jaccard.dist, "method") <- "jaccard"
    return(as.dist(jaccard.dist))
}

# Jaccard distance
hc.boot <- pvclust(df, method.hclust = 'ward.D2', method.dist = jaccard,
                   nboot = 1000, parallel = TRUE, iseed = 123456)

# Euclidean distance
#hc.boot <- pvclust(df, method.hclust = 'ward.D2', method.dist = 'euclidean',
#                  nboot = 1000, parallel = TRUE, iseed = 123456)

# Output files
#out_tree_file = paste(input_path, "/", input_name_noExt, ".tree", sep = "")
picture_type = "pdf"
out_picture_file = paste(input_folder, "/", "tree", ".", picture_type, sep = "")

# Open figure file handle
pdf(out_picture_file)

# Make tree
plot(hc.boot, cex = 0.3, cex.pv = 0.3)

# Add boxes around the significant clusters
#pvrect(hc.boot, alpha=0.95)

# Close picture file handle
dev.off()
