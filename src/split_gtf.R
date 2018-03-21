#!/usr/bin/env Rscript

library(GenomicFeatures)
library(rtracklayer)
library(systemPipeR)

###########
# GLOBALS #
###########

gtf_file <- snakemake@input[["gtf"]]
gff_file <- snakemake@input[["gff"]]
gene_output <- snakemake@output[["genic"]]
intergenic_output <- snakemake@output[["intergenic"]]
cpus <- snakemake@threads[[1]]
log_file <- snakemake@log[["log"]]

# dev
# gtf_file <- "output/010_ref/Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf"
# gff_file <- "data/ref/Araport11_GFF3_genes_transposons.201606.gff"
# cpus <- 8
# gtf_output <- paste0(
#     "output/010_ref/",
#     "Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf")

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# set up multiprocessing
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# prepare annotations
gff <- rtracklayer::import.gff(gff_file)
gtf <- rtracklayer::import.gff(gtf_file)

full_txdb <- GenomicFeatures::makeTxDbFromGRanges(
    gff[seqnames(gff) %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")]
)
txdb <- GenomicFeatures::makeTxDbFromGRanges(gtf)

# features to shuffle
all_genes <- genes(txdb)
all_intergenic <- genFeatures(full_txdb, "intergenic", reduce_ranges = FALSE)

# export features
export(all_genes, gene_output, "bed")
export(all_intergenic[[1]], intergenic_output, "bed")

# write session info
sessionInfo()
