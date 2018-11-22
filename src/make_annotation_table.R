#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(rtracklayer)

###########
# GLOBALS #
###########

gff_file <- snakemake@input[["gff"]]
annotation_file <- snakemake@output[["annotation"]]

# dev
# gff_file <- "data/ref/Araport11_GFF3_genes_transposons.201606.gff"


########
# MAIN #
########

# read gff
gff_genes <- import.gff3(gff_file, feature.type = "gene")

# parse annotations
gene_metadata <- as.data.table(mcols(gff_genes))
ara11_annotation <- gene_metadata[, .(ID, symbol, full_name, Note)]
setkey(ara11_annotation, ID)

# write 
fwrite(ara11_annotation, annotation_file)

# log
sessionInfo()

