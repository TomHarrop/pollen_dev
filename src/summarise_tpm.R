#!/usr/bin/env Rscript

log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

tpm_file <- snakemake@input[["tpm"]]
tpm_summary_file <- snakemake@output[["tpm_summary"]]
tpm_summary_wide_file <- snakemake@output[["tpm_summary_wide"]]
tpm_wide_file <- snakemake@output[["tpm_wide"]]

stage_order <- c("RUNM", "PUNM", "LBCP", "LTCP")

#dev 
# tpm_file <-  'output/080_filter-background/gene_calls.csv'

########
# MAIN #
########

# read the tpm file
tpm_long <- fread(tpm_file)

# remove first run
plants_to_keep <- paste0("p", 5:8)
tpm_subset <- tpm_long[plant %in% plants_to_keep]

# replace plant names with replicate names
rep_names <- c(p5 = "rep1", p6 = "rep2", p7 = "rep3", p8 = "rep4")
tpm_subset[, plant := plyr::revalue(plant, rep_names)]
tpm_subset[, stage := factor(stage, levels = stage_order)]

# make a wide data.table
tpm_wide <- dcast(tpm_subset, id ~ stage + plant, value.var = "tpm")

# remove undetected genes
tpm_summary_wide <- copy(tpm_subset)
tpm_summary_wide[detected_stage == FALSE, tpm := NA]

# summarise the tpm
tpm_summary <- tpm_summary_wide[, .(tpm_mean = mean(tpm),
             tpm_sd = sd(tpm),
             n = length(tpm[!is.na(tpm)])),
           by = .(stage, id)]
tpm_summary[, tpm_sem := tpm_sd / sqrt(n)]

col_order <- c("id",
               "tpm_mean_RUNM", 
               "tpm_sem_RUNM",
               "tpm_mean_PUNM",
               "tpm_sem_PUNM", 
               "tpm_mean_LBCP",
               "tpm_sem_LBCP", 
               "tpm_mean_LTCP", 
               "tpm_sem_LTCP")
tpm_summary_wide <- dcast(tpm_summary,
                          id ~ stage,
                          value.var = c("tpm_mean", "tpm_sem"))
setcolorder(tpm_summary_wide, col_order)

# write output
fwrite(tpm_summary_wide, tpm_summary_wide_file)
fwrite(tpm_summary, tpm_summary_file)
fwrite(tpm_wide, tpm_wide_file)

# record log
sessionInfo()
