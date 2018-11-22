#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)
library(DESeq2)
library(Mfuzz)
library(viridis)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds"]]
annotation_file <- snakemake@input[["annotation"]]

threads <- snakemake@threads[[1]]
alpha <- snakemake@params[["alpha"]]
seed <- snakemake@params[["seed"]]

cluster_plot_file <- snakemake@output[["cluster_plot"]]
annotated_clusters_file <- snakemake@output[["annotated_clusters"]]

stage_order <- c("RUNM", "PUNM", "LBCP", "LTCP")

# dev
# dds_file <- "output/090_deseq/dds.Rds"
# alpha <- 0.1
# threads <- 8
# seed <- 1
# annotation_file <- "output/010_ref/araport_annotation.csv"

########
# MAIN #
########

BiocParallel::register(BiocParallel::MulticoreParam(threads))

# read data
dds <- readRDS(dds_file)

# remove bad samples
plants_to_keep <- paste0("p", 5:8)
dds_subset <- dds[, dds$plant %in% plants_to_keep]
dds_subset$plant <- droplevels(dds_subset$plant)

# run a likelihood ratio test
dds_subset <- DESeq(dds_subset,
                    test = "LRT",
                    reduced = ~ plant,
                    parallel = TRUE)

# extract time-dependent genes
lrt_res <- results(dds_subset, alpha = alpha)
lrt_sig <- rownames(subset(lrt_res, padj < alpha))

# transform expression values
dds_sig <- dds_subset[lrt_sig, ]
vst <- varianceStabilizingTransformation(dds_sig, blind = FALSE)

# calculate mean VST
vst_long <- data.table(melt(assay(vst)))
vst_long[, c("stage", "plant") := tstrsplit(Var2, "_")]
vst_means <- vst_long[, .(mean_vst = mean(value)), by = .(Var1, stage)]
vst_means[, stage := factor(stage, levels = stage_order)]
expression_matrix <- as.matrix(data.frame(dcast(vst_means, Var1 ~ stage),
                                          row.names = "Var1"))

# set up mfuzz object
pheno_data <- data.frame(row.names = colnames(expression_matrix),
                         stage = factor(colnames(expression_matrix),
                                        levels = stage_order))

vg <- ExpressionSet(assayData = expression_matrix,
                    phenoData = new('AnnotatedDataFrame', data = pheno_data))

# standardise
vg_s <- standardise(vg)

# optimise paramaters
# m <- mestimate(vg_s) # too high
m <- 2.3
Dmin(vg_s, m, crange = seq(4, 10, 1), repeats = 3)

# run the clustering
message(paste("Clustering with seed", seed))
set.seed(seed)
c1 <- mfuzz(vg_s, c = 6, m = m)
clusters <- acore(vg_s, c1, min.acore = 0.7)

# quick check
mfuzz.plot(vg_s, c1,
           mfrow = c(2, 5),
           min.mem = 0.7,
           time.labels = levels(vg$stage),
           new.window = FALSE)
sapply(clusters, dim)

# annotate clusters
ara11 <- fread(annotation_file)
annotated_cluster_list <- lapply(clusters, function(x) 
    merge(x,
          ara11,
          all.x = TRUE,
          all.y = FALSE,
          by.x = "NAME",
          by.y = "ID")
)
annotated_clusters <- rbindlist(annotated_cluster_list, idcol = "cluster")

# combine data for plotting
cluster_membership <- rbindlist(clusters, idcol = "cluster")

cluster_expr_long <- data.table(exprs(vg_s), keep.rownames = TRUE)
setnames(cluster_expr_long, "rn", "NAME")
cluster_expr <- melt(cluster_expr_long,
                     id.vars = "NAME",
                     variable.name = "stage",
                     value.name = "scaled_vst")

cluster_pd <- merge(cluster_membership,
                    cluster_expr,
                    by = "NAME",
                    all.x = TRUE,
                    all.y = FALSE)

# re-order based on expression in RUNM
meanmean_expression <- cluster_pd[,
                                  .(stage_mean = mean(scaled_vst)),
                                  by = .(stage, cluster)]
meanmean_expression[, int_mean := round(stage_mean, 1)]
clust_ord <- dcast(meanmean_expression,
                   cluster ~ stage,
                   value.var = "int_mean")
setorder(clust_ord, -RUNM, -PUNM, -LBCP, -LTCP)
cluster_order <- 1:clust_ord[, length(unique(cluster))]
names(cluster_order) <- clust_ord[, unique(cluster)]

# RENAME THE CLUSTERS!!!
cluster_pd[, cluster := factor(plyr::revalue(factor(cluster),
                                      cluster_order),
                               levels = cluster_order)]
annotated_clusters[, cluster := factor(plyr::revalue(factor(cluster),
                                                     cluster_order),
                                       levels = cluster_order)]
setorder(annotated_clusters, cluster)

# plot the closes to core last?
setorder(cluster_pd, `MEM.SHIP`)
gene_order <- cluster_pd[, unique(NAME)]
cluster_pd[, NAME := factor(NAME, levels = gene_order)]

# label the clusters with the number of genes
cluster_pd[, genes_per_cluster := 
               length(unique(as.character(NAME))), by = cluster]
cluster_pd[, cluster_label := 
               paste0("Cluster ", cluster, " (", genes_per_cluster, " genes)")]
setorder(cluster_pd, cluster)
cluster_pd[, cluster_label :=
               factor(cluster_label, levels = unique(cluster_label))]

# plot clusters
gp <- ggplot(cluster_pd, aes(x = stage,
                       y = scaled_vst,
                       colour = `MEM.SHIP`,
                       group = NAME)) +
    theme_minimal(base_size = 8) +
    xlab(NULL) + ylab("Scaled, mapped reads") +
    facet_wrap(~ cluster_label) +
    scale_colour_viridis(guide = FALSE) +
    geom_path()

# write output
ggsave(cluster_plot_file, gp, width = 150, height = 125, units = "mm")
fwrite(annotated_clusters, annotated_clusters_file)

sessionInfo()
