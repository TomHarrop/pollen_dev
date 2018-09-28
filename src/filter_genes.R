#!/usr/bin/env Rscript

log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########

tpm_file <- snakemake@input[["tpm"]]
intergenic_tpm_file <- snakemake@input[["intergenic_tpm"]]
detected_genes_file <- snakemake@output[["detected_genes"]]
gene_call_file <- snakemake@output[["gene_calls"]]
gp1_file <- snakemake@output[["freqpoly"]]
gp2_file <- snakemake@output[["violin"]]

# dev
# tpm_file <- "output/070_tpm/tpm.csv"
# intergenic_tpm_file <- "output/060_cutoffs/intergenic_tpm.csv"
# detected_genes_file <- "test1.csv"
# gene_call_file <- "test2.csv"
# gp1_file <- "test1.pdf"
# gp2_file <- "test2.pdf"



########
# MAIN #
########

# combine data
intergenic_tpm <- fread(intergenic_tpm_file)
intergenic_tpm[, id := paste0(id, "_intergenic"), by = id]
tpm_data <- rbindlist(list("Genes" = fread(tpm_file),
                           "Intergenic" = intergenic_tpm),
                      idcol = "type",
                      use.names = TRUE)

# recalculate tpm
tpm_data[, c("T.g", "T.sum", "tpm") := NULL]
tpm_data[, T.g := as.numeric(counts) * as.numeric(rl) / as.numeric(feature_length),
         by = .(id, sample, type)]
tpm_data[, T.sum := sum(T.g), by = sample]
tpm_data[, tpm := (as.numeric(counts) * as.numeric(rl) * 1e6) / (feature_length * T.sum)]

# calculate cutoffs
tpm_data[, c("stage", "plant") := tstrsplit(sample, "_"), by = sample]
cutoffs <- tpm_data[, .(q95 = quantile(tpm, 0.95)),
                    by = .(sample, type)][type == "Intergenic", .(sample, q95)]
cutoffs[, c("stage", "plant") := tstrsplit(sample, "_"), by = sample]

# rearrange plot
tpm_data[, type := factor(type, levels = c("Intergenic", "Genes"))]
stage_order <- c("RUNM", "PUNM", "LBCP", "LTCP")
tpm_data[, stage := factor(stage, levels = stage_order)]
tpm_data[, plant := factor(plant, levels = sort(unique(plant)))]
cutoffs[, stage := factor(stage, levels = stage_order)]
cutoffs[, plant := factor(plant, levels = sort(unique(plant)))]

# filter
tpm_filter <- merge(tpm_data, cutoffs)
tpm_filter[type == "Genes", detected_lib := tpm > q95]
tpm_filter[type == "Genes",
           detected_stage := sum(detected_lib, na.rm = TRUE) >= 4,
           by = .(id, stage)]
detected_genes <- tpm_filter[detected_stage == TRUE, unique(id)]

# make plots
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
nudge <- 1
gp1 <- ggplot(tpm_data, aes(x = tpm + nudge, colour = type)) +
    scale_colour_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
    xlab(paste0("TPM + ", nudge)) + ylab("Count") +
    ggtitle("Frequency polygons (distributions), showing all genes") +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(stage ~ plant) +
    geom_vline(data = cutoffs, aes(xintercept = q95 + nudge),
               linetype = 2, colour = Set1[1]) +
    geom_freqpoly(alpha = 0.5, bins = 315)

gp2 <- ggplot(tpm_data, aes(y = tpm, x = type, fill = type)) +
    scale_fill_brewer(palette = "Set1", guide = FALSE) +
    xlab(NULL) + ylab("TPM") +
    ggtitle("Violin plots (distributions), only showing genes with TPM > 0") +
    scale_y_log10() +
    facet_grid(stage ~ plant) +
    geom_hline(data = cutoffs, aes(yintercept = q95),
               linetype = 2, colour = Set1[1]) +
    geom_violin(alpha = 0.5)

# write output
fwrite(data.table(detected_genes),
       detected_genes_file,
       col.names = FALSE)

fwrite(tpm_filter[type == "Genes", .(sample,
                              stage,
                              plant,
                              id,
                              counts,
                              tpm,
                              detected_lib,
                              detected_stage)],
       gene_call_file)
ggsave(gp1_file, gp1, width = 10, height = 7.5, units = "in")
ggsave(gp2_file, gp2, width = 10, height = 7.5, units = "in")

# write session info
sessionInfo()