#!/usr/bin/env Rscript

log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(DESeq2)

###########
# GLOBALS #
###########

gene_call_file <- snakemake@input[["gene_calls"]]
detected_gene_file <- snakemake@input[["detected_genes"]]
dds_file <- snakemake@output[["dds"]]

# dev
# gene_call_file <- 'output/080_filter-background/gene_calls.csv'
# detected_gene_file <- 'output/080_filter-background/detected_genes.csv'

########
# MAIN #
########

# read data
detected_genes <- fread(detected_gene_file, header = FALSE)[, unique(V1)]
gene_calls <- fread(gene_call_file)
det_calls <- gene_calls[id %in% detected_genes]

# generate count matrix
count_dt <- dcast(det_calls, id ~ sample, value.var = "counts")
count_matrix <- as.matrix(data.frame(count_dt, row.names = "id"))

# generate col data
stage_order <- c("UNM", "PUNM", "BCP", "TCP")
det_calls[, stage := factor(stage, levels = stage_order)]
det_calls[, plant := factor(plant, levels = sort(unique(plant)))]
cd <- data.frame(unique(det_calls[, .(sample, stage, plant)]),
                 row.names = "sample")


# generate deseq object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = cd,
                              design = ~ plant + stage)

# write output
saveRDS(dds, dds_file)

# write session info
sessionInfo()

