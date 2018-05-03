library(data.table)
library(ggplot2)

#############
# FUNCTIONS #
#############

ParseMu <- function(star_log_file){
    my_lines <- readLines(star_log_file)
    my_line <- grep("Average mapped length", my_lines, value = TRUE)
    as.numeric(gsub("[^[:digit:]\\.]+", "", my_line))
}

###########
# GLOBALS #
###########

count_files <- snakemake@input[["bg_counts"]]
star_dir <- snakemake@params[["star_dir"]]
tpm_output <- snakemake@output[["intergenic_tpm"]]
log_file <- snakemake@log[["log"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# get mu
star_log_files <- list.files(star_dir,
                             pattern = "Log.final.out",
                             full.names = TRUE)
names(star_log_files) <- gsub("^([^\\.]+).*", "\\1", basename(star_log_files))
mu_values <- lapply(star_log_files, ParseMu)
mu_table <- data.table(rl = mu_values,
           sample = names(mu_values))

# get counts
count_list <- lapply(count_files, fread)
counts <- rbindlist(count_list)

# merge counts and mu
tpm <- merge(counts, mu_table, by = "sample", all = TRUE)

tpm[, T.g := as.numeric(counts) * as.numeric(rl) / as.numeric(feature_length),
    by = .(id, sample)]
tpm[, T.sum := sum(T.g), by = sample]
tpm[, tpm := (as.numeric(counts) * as.numeric(rl) * 1e6) / (feature_length * T.sum)]

# write output
fwrite(tpm, tpm_output)

# write session info
sessionInfo()

 