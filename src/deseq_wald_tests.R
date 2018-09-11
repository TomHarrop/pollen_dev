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

threads <- snakemake@threads[[1]]

alpha <- snakemake@params[["alpha"]]
lfcThreshold <- snakemake@params[["lfc_threshold"]]

# dev
# dds_file <- "output/090_deseq/dds.Rds"
# alpha <- 0.1
# lfcThreshold <- log(1.5, 2)

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
               c("UNM", "PUNM"), "UNM_PUNM", "BCP_TCP"),
    levels = c("UNM_PUNM", "BCP_TCP"))
design(dds_group) <-  ~ plant + group

# run the wald tests
dds_subset <- DESeq(dds_subset, parallel = TRUE)
dds_group <- DESeq(dds_group, parallel = TRUE)

# extract group results
group_res <- results(dds_group,
                     name = "group_BCP_TCP_vs_UNM_PUNM",
                     lfcThreshold = lfcThreshold,
                     alpha = alpha,
                     tidy = TRUE)
group_dt <- data.table(group_res)
setnames(group_dt, "row", "id")
group_dt[, wald_test := "group_BCP_TCP_vs_UNM_PUNM"]

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

# write data
fwrite(group_dt, group_dt_file)
fwrite(pairwise_results, pairwise_results_file)


