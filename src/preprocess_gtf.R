#!/usr/bin/env Rscript

library(rtracklayer)

###########
# GLOBALS #
###########

gtf_file <- snakemake@input[["gtf"]]
gff_file <- snakemake@input[["gff"]]
gtf_output <- snakemake@output[["gtf"]]
cpus <- snakemake@threads[[1]]
log_file <- snakemake@log[["log"]]


# gtf_file <- "data/ref/Araport11_GFF3_genes_transposons.201606.gtf"
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
gtf <- rtracklayer::import.gff(gtf_file, format = "gtf")

# found the nuclear rRNA genes
ribo_gff <- subset(gff, type == "rRNA" & !seqnames %in% c("ChrC", "ChrM"))
ribo_parents <- unique(unlist(ribo_gff$Parent))

# get all nuclear annotations
nuclear_gtf <- subset(gtf, !seqnames %in% c("ChrC", "ChrM")) 

# eliminate nuclear rRNA genes
chuck <- apply(sapply(ribo_parents, function(x)
    grepl(x, nuclear_gtf$transcript_id)), 1, any)
nuc_no_rRNA <- nuclear_gtf[!chuck]

# write output
export.gff(nuc_no_rRNA, gtf_output, format = "gtf")

# write session info
sessionInfo()
