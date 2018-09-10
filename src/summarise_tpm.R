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
tpm_wide_file <- snakemake@output[["tpm_wide"]]

########
# MAIN #
########

# read the tpm file
tpm_long <- fread(tpm_file)
tpm_long[, c("stage", "plant") := tstrsplit(sample, "_")]

# remove first run
plants_to_keep <- paste0("p", 5:8)
tpm_subset <- tpm_long[plant %in% plants_to_keep]

# make a wide data.table
tpm_wide <- dcast(tpm_subset, id ~ sample, value.var = "tpm")

# summarise the tpm
tpm_summary <- tpm_subset[, .(tpm_mean = mean(tpm),
             tpm_sd = sd(tpm),
             n = length(tpm)),
           by = .(stage, id)]
tpm_summary[, tpm_sem := tpm_sd / sqrt(n)]

# write output
fwrite(tpm_summary, tpm_summary_file)
fwrite(tpm_wide, tpm_wide_file)

# record log
sessionInfo()
