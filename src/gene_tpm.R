library(data.table)
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)

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

gtf_file <- snakemake@input[["gtf"]]
count_files <- snakemake@input[["count_files"]]
star_dir <- snakemake@params[["star_dir"]]
tpm_output <- snakemake@output[["tpm"]]
log_file <- snakemake@log[["log"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# get feature lengths
gtf <- rtracklayer::import.gff(gtf_file)
grl <- reduce(split(gtf, elementMetadata(gtf)$gene_id))
gtf.reduced <- unlist(grl, use.names = FALSE)
elementMetadata(gtf.reduced)$gene_id <- rep(names(grl), elementNROWS(grl))
elementMetadata(gtf.reduced)$widths <- width(gtf.reduced)
feature_lengths <- data.table(data.frame(elementMetadata(gtf.reduced)))[
    ,.(feature_length = sum(widths)), by = gene_id]
setnames(feature_lengths, "gene_id", "id")

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
names(count_list) <- gsub("^([^\\.]+).*", "\\1", basename(count_files))
counts <- rbindlist(count_list, idcol = "sample")[
    , .(id = V1, counts = V2, sample)][!grepl("^N_", id)]

# merge all
tpm <- merge(merge(counts, feature_lengths),
      mu_table, by = "sample")

# calculate tpm
tpm[, T.g := as.numeric(counts) * as.numeric(rl) / as.numeric(feature_length),
    by = .(id, sample)]
tpm[, T.sum := sum(T.g), by = sample]
tpm[, tpm := (as.numeric(counts) * as.numeric(rl) * 1e6) /
        (feature_length * T.sum)]

# write output
fwrite(tpm, tpm_output)

# write session info
sessionInfo()
