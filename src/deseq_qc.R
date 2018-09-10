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
log_file <- snakemake@log[["log"]]

# dev
# dds_file <- "output/090_deseq/dds.Rds"

########
# MAIN #
########

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

# plot the pca
pc_dt[, c("stage", "plant") := tstrsplit(sample_id, "_")]
stage_order <- c("UNM", "PUNM", "BCP", "TCP")
pc_dt[, stage := factor(stage, levels = stage_order)]
pc_dt[, plant := factor(plant, levels = sort(unique(plant)))]

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
    facet_wrap(~ facet_label, nrow = 2) +
    scale_colour_brewer(palette = "Set1") +
    geom_point(size = 2,
               alpha = 0.7,
               shape = 16,
               position = position_jitter(width = 0.2))

ggsave("pca.pdf", gp, width = 10, height = 7.5, units = "in")

# plot a heatmap
colours <- colorRampPalette(rev(RColorBrewer::brewer.pal(6, "Blues")))(255)
sample_dists <- dist(t(assay(vst)), method = "minkowski")
cairo_pdf("heatmap.pdf", width = 10, height = 7.5)
pheatmap::pheatmap(as.matrix(sample_dists),
                   clustering_distance_cols = sample_dists,
                   clustering_distance_rows = sample_dists,
                   col = colours)
dev.off()

# anything significant?
subset(results(dds_lrt), padj < 0.05)
res_lrt <- results(dds_lrt)
rownames(head(res_lrt[order(res_lrt$padj), ]))

PlotMyGene <- function(goi) {
    my_counts <- counts(dds_lrt, normalized = TRUE)[goi,]
    my_count_dt <- data.table(sample_id = names(my_counts),
                              normalized_count = my_counts)
    stage_order <- c("UNM", "PUNM", "BCP", "TCP")
    my_count_dt[, c("stage", "plant") := tstrsplit(sample_id, "_")]
    my_count_dt[, stage := factor(stage, levels = stage_order)]
    my_count_dt[, plant := factor(plant, levels = sort(unique(plant)))]
    ggplot(my_count_dt, aes(x = stage,
                            y = normalized_count,
                            colour = plant,
                            group = 1)) +
        scale_colour_brewer(palette = "Set1") +
        scale_y_log10() +
        geom_smooth(colour = alpha("black", 0.5),
                    se = FALSE,
                    span = 0.5) +
        geom_point(size = 2,
                   alpha = 0.7,
                   shape = 16,
                   position = position_jitter(width = 0.2))
}

x <- counts(dds_lrt, normalized = TRUE)["AT5G09330",]
fwrite(data.table(sample = names(x), normalised_counts = x), "AT5G09330.csv")
my_plot <- PlotMyGene("AT5G09330")
ggsave("AT5G09330.pdf", my_plot, width = 5, height = 3.5, units = "in")

# try a wald test
dds_wald <- DESeq(dds,
                  parallel = TRUE)
res_wald1 <- results(dds_wald,
                     name = "stage_TCP_vs_UNM",
                     lfcThreshold = log(2, 2),
                     alpha = 0.05)
subset(res_wald1, padj < 0.05)

res_wald2 <- results(dds_wald,
                     name = "stage_PUNM_vs_UNM",
                     lfcThreshold = log(2, 2),
                     alpha = 0.05)
subset(res_wald2, padj < 0.05)
subset(res_wald2[order(res_wald2$padj),], baseMean > 500)


my_plot <- PlotMyGene("AT1G23570")
ggsave("gene.pdf", device = cairo_pdf(), my_plot, width = 10, height = 7.5, units = "in")


plotCounts(dds_wald, "AT3G29078", intgroup = c("stage"), returnData = TRUE)

# specific gene counts
bpcs <- c("AT2G01930", "AT1G14685", "AT1G68120", "AT2G21240", "AT4G38910", 
          "AT5G42520", "AT2G35550")
wcounts <- data.table(data.frame(counts(dds_wald, normalized = TRUE)),
                      keep.rownames = TRUE)
pd_wide <- wcounts[rn %in% bpcs]

pd <- melt(pd_wide, id.vars = "rn", variable.name = "sample", value.name = "norm_count")
pd[, c("stage", "plant") := tstrsplit(sample, "_")]
stage_order <- c("UNM", "PUNM", "BCP", "TCP")
pd[, stage := factor(stage, levels = stage_order)]
gp <- ggplot(pd, aes(y = norm_count, x = stage, group = 1)) +
    xlab(NULL) + ylab("Normalized read count") +
    scale_colour_brewer(palette = "Set1") +
    facet_wrap(~ rn, scales = "free_y") +
    geom_smooth(method = "glm", se = FALSE, colour = alpha("blue", 0.5)) +
    geom_point(shape = 16, position = position_jitter(width = 0.1),
               mapping = aes(colour = plant))
ggsave("BPCs.pdf", gp, width = 10, height = 7.5, units = "in")

op <- dcast(pd, rn ~ stage+plant, value.var = "norm_count")
setnames(op, "rn", "gene_id")
fwrite(op, "bpc.csv")

# why are we losing all these genes?
kept_counts <- counts(dds_filtered, normalized = TRUE)
discarded_genes <- rownames(dds_lrt[is.na(results(dds_lrt)$padj), ])
lost_counts <- counts(dds_lrt[discarded_genes, ], normalized = TRUE)
kept_lost_wide <- rbindlist(lapply(list(
    kept = kept_counts,
    lost = lost_counts),
    data.table, keep.rownames = TRUE), idcol = "type")
kept_lost <- melt(kept_lost_wide, id.vars = c("type", "rn"),
                  variable.name = "sample",
                  value.name = "norm_counts")
gene_var <- kept_lost[, .(gene_var = var(norm_counts)), by = .(type, rn)]

# because they have low counts
ggplot(kept_lost, aes(x = type, y = norm_counts + 0.5)) +
    scale_y_log10() +
    geom_violin() +
    facet_wrap(~ sample, scales = "free_y") +
    geom_hline(yintercept = 1)

# and because they have large variance
ggplot(gene_var, aes(x = type, y = gene_var)) +
    scale_y_log10() +
    geom_violin()

plotMA(dds_lrt)
plotDispEsts(dds_lrt, )
plotDispEsts(dds_filtered, )

# transform
vst <- varianceStabilizingTransformation(dds_filtered, blind = TRUE)

# select genes to assess
ntop <- Inf
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                   length(rv)))]

# principal components
pca_res <- prcomp(t(assay(vst)[select, ]), center = TRUE)
pv <- pca_res$sdev^2/sum(pca_res$sdev^2)

# percent var
pv_dt <- data.table(PC = colnames(pca_res$rotation),
                    pv = pv)
pv_dt[, percent_var := round(pv * 100, 1)]

# extract pc results
pca_dt <- data.table(pca_res$x, keep.rownames = TRUE)
cd <- data.table(as.data.frame(colData(dds)), keep.rownames = TRUE)
pca_wide <- merge(pca_dt, cd)

# generate plot data
pca_long <- melt(pca_wide,
                 id.vars = c("stage", "plant", "rn"),
                 variable.name = "PC",
                 value.name = "score")
pca_pd <- merge(pca_long, pv_dt)
pca_pd[, pc_label := paste0(PC, " (", percent_var, "%)")]
pca_pd[, pc_num := as.numeric(gsub("[^[:digit:]]+", "", PC))]
setorder(pca_pd, pc_num)
pca_pd[, pc_label := factor(pc_label, levels = unique(pc_label))]

# draw plot
ggplot(pca_pd, aes(x = stage, y = score, colour = plant)) +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~ pc_label) +
    geom_point(position = position_jitter(width = 0.1), shape = 16)

plotPCA(vst, intgroup = "stage", ntop = 3000)
