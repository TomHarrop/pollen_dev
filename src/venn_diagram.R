# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(readxl)
library(data.table)
library(VennDiagram)

###########
# GLOBALS #
###########

array_file <- snakemake@input[["array"]]
call_file <- snakemake@input[["calls"]]
tapetum_file <- snakemake@input[["tapetum"]]
venn_diagram_file <- snakemake@output[["venn_diagram"]]
array_comparison_file <- snakemake@output[["array_comparison"]]
det_in_file <- snakemake@output[["detected_genes_matrix"]]

########
# MAIN #
########

# read data
array_data_tb <- read_xls(array_file)
array_data <- as.data.table(array_data_tb)
call_data <- fread(call_file)
tapetum_data <- read_xls(tapetum_file,
         sheet = 2,
         skip = 5,
         col_names = c("id", "tap8-10"))


# list array genes
array_unm <- array_data[UNM > 0, unique(`Gene Name`)]
array_bcp <- array_data[BCP > 0, unique(`Gene Name`)]
array_tcp <- array_data[TCP > 0, unique(`Gene Name`)]

# list tapetum genes
tapetum_genes <- unique(tapetum_data$id)
tapetum_genes <- tapetum_genes[!is.na(tapetum_genes)]

# list rna seq genes
detected_stage <- dcast(unique(call_data, by = c("stage", "id", "detected_stage")),
                        id ~ stage,
                        value.var = "detected_stage")
pol_lbcp <- detected_stage[LBCP == TRUE, unique(id)]
pol_ltcp <- detected_stage[LTCP == TRUE, unique(id)]
pol_runm <- detected_stage[RUNM == TRUE, unique(id)]
pol_punm <- detected_stage[PUNM == TRUE, unique(id)]

# list special RNA seq combo
pol_runmpunm <- unique(union(pol_runm, pol_punm))

# make a list of all lists
intersect_list <- list(RUNM = pol_runm,
                       PUNM = pol_punm,
                       LBCP = pol_lbcp,
                       LTCP = pol_ltcp,
                       `Array\nUNM` = toupper(array_unm),
                       `Array\nBCP` = toupper(array_bcp),
                       `Array\nTCP` = toupper(array_tcp),
                       `RUNM &\nPUNM` = pol_runmpunm,
                       Tapetum = tapetum_genes)

# per gene "detected in" table
all_genes <- unique(unlist(intersect_list))
det_in_list <- lapply(all_genes, function(x)
    lapply(intersect_list, function(y)
        x %in% y))
names(det_in_list) <- all_genes
det_in <- rbindlist(det_in_list, idcol = "id")
names(det_in) <- sub("\n", "_", names(det_in))
names(det_in) <- sub(" ", "_", names(det_in))
fwrite(det_in, det_in_file)

# sapply(intersect_list, function(x) length(unique(x)))
# length(unique(intersect(pol_runmpunm, toupper(array_unm))))
# length(unique(intersect(pol_runmpunm, toupper(array_unm))))/length(unique(array_unm))
# 
# length(unique(intersect(pol_lbcp, toupper(array_bcp))))
# length(unique(intersect(pol_ltcp, toupper(array_tcp))))

# make a comparison table for array vs. rnaseq
array_comparison <- data.table(
    stage = c("RUNM and PUNM", "BCP", "TCP"),
    rna_seq = c(length(pol_runmpunm), length(pol_lbcp), length(pol_ltcp)),
    array = c(length(array_unm), length(array_bcp), length(array_tcp)),
    common = c(length(unique(intersect(pol_runmpunm, toupper(array_unm)))),
               length(unique(intersect(pol_lbcp, toupper(array_bcp)))),
               length(unique(intersect(pol_ltcp, toupper(array_tcp))))))
array_comparison[, percent_detected_in_rnaseq := round(common * 100 / array, 1)]
array_comparison[, number_unique_to_rnaseq := rna_seq - common]
fwrite(array_comparison, array_comparison_file)

# print the venn diagrams
Set1 <- RColorBrewer::brewer.pal(9, "Set1")

pdf(venn_diagram_file, width = 4, height = 4, pointsize = 8)

# rnaseq dataset venn diagram
vd1 <- venn.diagram(
    intersect_list[c(1, 4, 2, 3)],
    filename = NULL,
    fill = Set1[c(1:4)],
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 1,
    fontfamily = 'Helvetica',
    cat.fontfamily = 'Helvetica',
    alpha = 0.5,
    margin = 0)
grid.draw(vd1)

grid.newpage()
vd2 <- venn.diagram(
    intersect_list[c(8, 5)],
    filename = NULL,
    fill = Set1[c(1:2)],
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 1,
    fontfamily = 'Helvetica',
    cat.fontfamily = 'Helvetica',
    alpha = 0.5,
    margin = 0.1,
    cat.just = list(c(1, 0.5), c(0, 0.5)))
grid.draw(vd2)

grid.newpage()
vd3 <- venn.diagram(
    intersect_list[c(3, 4, 6, 7)],
    fill = Set1[c(1:4)],
    filename = NULL,
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 1,
    fontfamily = 'Helvetica',
    cat.fontfamily = 'Helvetica',
    alpha = 0.5,
    distance = 0.1)

grid.draw(vd3)

grid.newpage()

vd4 <- venn.diagram(
    intersect_list[c(1, 2, 9)],
    fill = Set1[c(1:3)],
    filename = NULL,
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 1,
    fontfamily = 'Helvetica',
    cat.fontfamily = 'Helvetica',
    alpha = 0.5,
    distance = 0.1)

grid.draw(vd4)

dev.off()


# log
sessionInfo()

# print an upset for pollen_dev data
# can't use upset, they don't obey graphic parameters
# upset(fromList(intersect_list[c(1:4)]),
#       sets = rev(names(intersect_list)[c(1:4)]),
#       nsets = 6,
#       keep.order = TRUE,
#       empty.intersections = "on",
#       sets.x.label = "Number of genes",
#       point.size = 1)
# 
# # rnaseq vs early array
# upset(fromList(intersect_list[c(5, 8)]),
#       sets = rev(names(intersect_list)[c(5, 8)]),
#       nsets = 6,
#       keep.order = TRUE,
#       empty.intersections = "on",
#       sets.x.label = "Number of genes",
#       point.size = 1)
# 
# # rnaseq ltcp / lbcp vs array
# upset(fromList(intersect_list[c(3,4,6,7)]),
#       sets = rev(names(intersect_list)[c(3,4,6,7)]),
#       nsets = 6,
#       keep.order = TRUE,
#       empty.intersections = "on",
#       sets.x.label = "Number of genes")
# 
# 
# upset(fromList(intersect_list[c(1:3)]),
#       sets = rev(names(intersect_list)[c(1:3)]),
#       nsets = 6,
#       keep.order = TRUE,
#       empty.intersections = "on",
#       sets.x.label = "Number of genes")
