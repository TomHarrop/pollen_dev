#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(DESeq2)

#############
# FUNCTIONS #
#############

GetStageContrastResults <- function(contrast, name, dds) {
    my_contrast <- c(name, contrast[[2]], contrast[[1]])
    my_results <- data.table(results(dds,
                                     contrast = my_contrast,
                                     lfcThreshold = lfcThreshold,
                                     alpha = alpha,
                                     tidy = TRUE))
    setnames(my_results, "row", "id")
    my_results[, wald_test := paste(contrast, collapse = "_vs_")]
    return(my_results)
}

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds"]]

group_dt_file <- snakemake@output[["group_test"]]
pairwise_results_file <- snakemake@output[["stage_tests"]]
de_matrix_file <- snakemake@output[["de_matrix"]]

threads <- snakemake@threads[[1]]

alpha <- snakemake@params[["alpha"]]
lfcThreshold <- snakemake@params[["lfc_threshold"]]

stage_order <- c("RUNM", "PUNM", "LBCP", "LTCP")

# dev
# dds_file <- "output/090_deseq/dds.Rds"
# alpha <- 0.1
# lfcThreshold <- log(1.5, 2)
# threads <- 8

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

# group for (UNM + pUNM) vs (BCP + TCP) comparison
dds_group <- copy(dds_subset)
dds_group$group <- factor(
    ifelse(as.character(dds_group$stage) %in%
               c("RUNM", "PUNM"), "RUNM_PUNM", "LBCP_LTCP"),
    levels = c("RUNM_PUNM", "LBCP_LTCP"))
design(dds_group) <-  ~ plant + group

# run the wald tests
dds_subset <- DESeq(dds_subset, parallel = TRUE)
dds_group <- DESeq(dds_group, parallel = TRUE)

# extract group results
group_res <- results(dds_group,
                     name = "group_LBCP_LTCP_vs_RUNM_PUNM",
                     lfcThreshold = lfcThreshold,
                     alpha = alpha,
                     tidy = TRUE)
group_dt <- data.table(group_res)
setnames(group_dt, "row", "id")
group_dt[, wald_test := "group_LBCP_LTCP_vs_RUNM_PUNM"]

# extract pairwise stage results
combos <- combn(unique(as.character(dds_subset$stage)),
                m = 2, paste, collapse = "_")
contrasts <- lapply(combos, function(x) unlist(strsplit(x, split = "_")))
names(contrasts) <- combos

pairwise_results_list <- lapply(contrasts,
                                GetStageContrastResults,
                                name = "stage",
                                dds = dds_subset)
pairwise_results <- rbindlist(pairwise_results_list)

# make pairwise comparison table
n_de_table <- pairwise_results[padj < alpha,
                               .(n_de = length(unique(id))),
                               by = wald_test]
n_de_table[, c("stage_1", "stage_2") := tstrsplit(wald_test, "_vs_")]
n_de_table[, wald_test := NULL]

all_comparisons <- rbind(rbind(data.table(
    n_de = NA,
    stage_1 = stage_order,
    stage_2 = stage_order),
    n_de_table[, .(n_de, stage_1 = stage_2, stage_2 = stage_1)]),
    n_de_table)

all_comparisons[, stage_1 := factor(stage_1, levels = stage_order)]
all_comparisons[, stage_2 := factor(stage_2, levels = stage_order)]

de_matrix <- as.matrix(
    data.frame(
        dcast(all_comparisons, stage_1 ~ stage_2, value.var = "n_de"),
        row.names = "stage_1"))

# write data
fwrite(group_dt, group_dt_file)
fwrite(pairwise_results, pairwise_results_file)
write.csv(de_matrix, de_matrix_file, quote = FALSE)

sessionInfo()
