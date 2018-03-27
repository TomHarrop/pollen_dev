#!/usr/bin/env Rscript

library(BiocParallel)
library(data.table)
library(GenomicAlignments)
library(rtracklayer)

###########
# GLOBALS #
###########

bamfile <- snakemake@input[["bam"]]
regions_file <- snakemake@input[["regions"]]
counts_file <- snakemake@output[["counts"]]
log_file <- snakemake@log[["log"]]

# dev
# bamfile <- "output/030_star-pass2/TCP_p2.Aligned.sortedByCoord.out.bam"
# regions_file <- "test/shuffled_valr.gff3"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# load regions
features <- import.gff3(regions_file)

# load bamfile
bam <- BamFile(bamfile)

# count
counts <- summarizeOverlaps(features = features,
                            reads = bam,
                            mode = "Union",
                            singleEnd = FALSE,
                            ignore.strand = TRUE,
                            fragments = TRUE)

# generate output
output_data <- data.table(
    sample = gsub("^([^\\.]+).*", "\\1", basename(bamfile)),
    counts = assay(counts)[, 1],
    id = rowData(counts)$ID,
    feature_length = width(counts)
)

# write output
fwrite(output_data, counts_file)

# write session info
sessionInfo()

