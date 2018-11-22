log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

stage_tests_file <- snakemake@input[["stage_tests"]]
cycle_genes_file <- snakemake@input[["cycle_genes"]]
annot_file <- snakemake@input[["annot"]]
cycle_wald_file <- snakemake@output[["cycle_wald"]]

# dev
# stage_tests_file = 'output/090_deseq/wald_stage.csv'
# cycle_genes_file = 'test.txt'
# cycle_wald_file = 'output/090_deseq/runm_punm_cell_cyle.csv'
# annot_file = 'output/010_ref/araport_annotation.csv'

########
# MAIN #
########

# read data
stage_tests <- fread(stage_tests_file)
cycle_genes <- fread(cycle_genes_file, header = FALSE, col.names = "id")
annot <- fread(annot_file)

# subset
cycle_results <- merge(stage_tests[wald_test == "PUNM_vs_RUNM"],
      cycle_genes)
annot_cycle_results <- merge(cycle_results, annot, by.x = "id", by.y = "ID")

# output
fwrite(annot_cycle_results, cycle_wald_file)

# log
sessionInfo()

