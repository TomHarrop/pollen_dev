#!/usr/bin/env Rscript

library(GenomicFeatures)
library(rtracklayer)
library(systemPipeR)
library(valr)

###########
# GLOBALS #
###########

gtf_file <- snakemake@input[["gtf"]]
gff_file <- snakemake@input[["gff"]]
seqlengths_file <- snakemake@input[["seqlengths"]]
shuffled_output <- snakemake@output[["shuffled"]]
cpus <- snakemake@threads[[1]]
log_file <- snakemake@log[["log"]]

# dev
# gtf_file <- "output/010_ref/Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf"
# gff_file <- "data/ref/Araport11_GFF3_genes_transposons.201606.gff"
# seqlengths_file <- "output/040_shuffle/seqlen.txt"
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
full_txdb <- GenomicFeatures::makeTxDbFromGRanges(
    gff[seqnames(gff) %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")]
)

gtf <- rtracklayer::import.gff(gtf_file)
txdb <- GenomicFeatures::makeTxDbFromGRanges(gtf)

# features to shuffle
all_genes <- genes(txdb)
width_cutoff <- quantile(width(all_genes), 0.99)
short_genes <- all_genes[width(all_genes) < width_cutoff]
all_intergenic <- genFeatures(full_txdb, "intergenic", reduce_ranges = FALSE)

# shuffle with valr
genome <- read_genome(seqlengths_file)
short_genes_tbl <- as.tbl_interval(short_genes)
short_genes_tbl$name <- names(short_genes)
shuffled_tbl <- bed_shuffle(short_genes_tbl,
                            genome = genome,
                            incl = as.tbl_interval(all_intergenic[[1]]),
                            within = TRUE,
                            seed = 7)

# setup for export
shuffled_gr <- makeGRangesFromDataFrame(shuffled_tbl,
                                        keep.extra.columns=FALSE,
                                        ignore.strand=TRUE,
                                        seqinfo=NULL,
                                        seqnames.field="chrom",
                                        start.field="start",
                                        end.field="end",
                                        strand.field="strand")
names(shuffled_gr) <- shuffled_tbl$name

# export features
export(shuffled_gr, shuffled_output, "gff3")

# write session info
sessionInfo()
