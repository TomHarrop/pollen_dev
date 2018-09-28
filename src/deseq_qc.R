#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(DESeq2)
library(ggplot2)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds"]]

pca_plot <- snakemake@output[["pca_plot"]]
distance_heatmap <- snakemake@output[["distance_heatmap"]]
pca_dt <- snakemake@output[["pca_dt"]]

threads <- snakemake@threads[[1]]

# dev
# dds_file <- "output/090_deseq/dds.Rds"

########
# MAIN #
########

BiocParallel::register(BiocParallel::MulticoreParam(threads))

# run the VST
dds <- readRDS(dds_file)

plants_to_keep <- paste0("p", 5:8)
dds_subset <- dds[, dds$plant %in% plants_to_keep]
dds_subset$plant <- droplevels(dds_subset$plant)

# extract genes with time-dependency
dds_lrt <- DESeq(dds_subset,
                 test = "LRT",
                 reduced = ~ plant,
                 parallel = TRUE)
keep <- !is.na(results(dds_lrt)$padj)
dds_filtered <- dds_lrt[keep, ]

# make a pca
vst <- varianceStabilizingTransformation(dds_filtered, blind = TRUE)
pc <- prcomp(t(assay(vst)), center = TRUE, scale = TRUE)
percent_exp <- pc$sdev^2/sum(pc$sdev^2) * 100

pc_dt <- data.table(pc$x, keep.rownames = TRUE)
setnames(pc_dt, "rn", "sample_id")

# replace plant names with replicate names
pc_dt[, c("stage", "plant") := tstrsplit(sample_id, "_")]
rep_names <- c(p5 = 1, p6 = 2, p7 = 3, p8 = 4)
pc_dt[, plant := plyr::revalue(plant, rep_names)]

# plot the pca
stage_order <- c("RUNM", "PUNM", "LBCP", "LTCP")
pc_dt[, stage := factor(stage, levels = stage_order)]
#pc_dt[, plant := factor(plant, levels = sort(unique(plant)))]

pc_long <- melt(pc_dt,
                id.vars = c("stage", "plant"),
                measure.vars = paste0("PC", 1:4),
                variable.name = "component")
lab_dt <- data.table(component = paste0("PC", 1:4),
                     percent_exp = percent_exp[1:4])
lab_dt[, facet_label := paste0(component,
                               " (",
                               round(percent_exp, 1),
                               "%)")]

pc_pd <- merge(pc_long, lab_dt)

gp <- ggplot(pc_pd, aes(x = stage, y = value, colour = plant)) +
    theme_minimal(base_size = 8) +
    xlab(NULL) + ylab("Score") +
    facet_wrap(~ facet_label, nrow = 2) +
    scale_colour_brewer(palette = "Set1",
                        guide = guide_legend(title = "Replicate")) +
    geom_point(size = 2,
               alpha = 0.7,
               shape = 16,
               position = position_jitter(width = 0.2))


# plot a heatmap
colours <- colorRampPalette(rev(RColorBrewer::brewer.pal(6, "Blues")))(255)
sample_dists <- dist(t(assay(vst)), method = "minkowski")


# write the output
fwrite(pc_dt, pca_dt)
ggsave(pca_plot, gp, width = 150, height = 125, units = "mm")
pheatmap::pheatmap(as.matrix(sample_dists),
                       clustering_distance_cols = sample_dists,
                       clustering_distance_rows = sample_dists,
                       col = colours,
                       filename = distance_heatmap)

# log
sessionInfo()
