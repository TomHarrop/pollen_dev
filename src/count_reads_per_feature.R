#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(data.table)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)
library(Rsamtools)
library(rtracklayer)
library(systemPipeR)

#############
# FUNCTIONS #
#############

FirstElement <- function(y){unlist(y)[1][1]}

ReadCountWrapper <- function(features,
                             bamfile_list) {
    my_counts <- GenomicAlignments::summarizeOverlaps(
        features = features,
        reads = bamfile_list,
        mode = "Union",
        singleEnd = TRUE,
        ignore.strand = TRUE,
        fragments = FALSE)
    my_rd <- rowData(my_counts)
    my_unlisted_rd_list <- lapply(seq_len(ncol(my_rd)), function(i)
        sapply(my_rd[,i], FirstElement))
    names(my_unlisted_rd_list) <- names(my_rd)
    my_unlisted_rd <- data.table(do.call(cbind, my_unlisted_rd_list))
    data.table(my_unlisted_rd,
               data.table(assay(my_counts)))
}

###########
# GLOBALS #
###########

bamfile_list <- snakemake@input[["bam_files"]]
gff_file <- snakemake@input[["gff"]]
plot1_file <- snakemake@output[["counts_plot"]]
plot2_file <- snakemake@output[["intron_exon_plot"]]
fc_file <- snakemake@output[["feature_counts"]]
cpus <- snakemake@threads[[1]]


#  cpus <- 6
# gff_file = "data/ref/Araport11_GFF3_genes_transposons.201606.gff"
# bamfile_list <- c("output/030_star-pass2/BCP_p1.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/BCP_p2.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/BCP_p3.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/BCP_p4.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/BCP_p5.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/BCP_p6.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/BCP_p7.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/BCP_p8.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p1.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p2.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p3.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p4.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p5.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p6.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p7.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/PUNM_p8.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p1.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p2.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p3.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p4.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p5.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p6.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p7.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/TCP_p8.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p1.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p2.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p3.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p4.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p5.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p6.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p7.Aligned.sortedByCoord.out.bam", 
#                   "output/030_star-pass2/UNM_p8.Aligned.sortedByCoord.out.bam")
# bamfile_list <- c("output/030_star-pass2/BCP_p5.Aligned.sortedByCoord.out.bam",
#                   "output/030_star-pass2/BCP_p6.Aligned.sortedByCoord.out.bam")
                  
#  plot1_file <- "counts_per_category.pdf"
#  plot2_file <- "intron_exon_counts.pdf"

########
# MAIN #
########

# set up multiprocessing
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# load bamfiles
names(bamfile_list) <- gsub("\\..*$", "", basename(bamfile_list))
bamfiles <- BamFileList(bamfile_list)

# prepare annotations
gff <- rtracklayer::import.gff(gff_file)

# found the nuclear rRNA genes
ribo_gff <- subset(gff, type == "rRNA" & !seqnames %in% c("ChrC", "ChrM"))
ribo_parents <- unique(unlist(ribo_gff$Parent))

# get all nuclear annotations
nuclear_gff <- subset(gff, !seqnames %in% c("ChrC", "ChrM")) 
nuclear_txdb_with_rrna <- GenomicFeatures::makeTxDbFromGRanges(nuclear_gff)

# remove rRNA genes from nuclear annotations
chuck <- apply(sapply(ribo_parents, function(x)
    grepl(x, nuclear_gff$ID)), 1, any)
nuc_no_rRNA <- nuclear_gff[!chuck]

# generate nuclear-only non-rRNA txdb
nuclear_txdb <- GenomicFeatures::makeTxDbFromGRanges(nuc_no_rRNA)
nuc_exons <- GenomicFeatures::exonicParts(nuclear_txdb,
                                          linked.to.single.gene.only = TRUE)
nuc_introns <- GenomicFeatures::intronicParts(nuclear_txdb,
                                              linked.to.single.gene.only = TRUE)

# mito regions
mito_gff <- subset(gff, seqnames == "ChrM")
mito_txdb <- GenomicFeatures::makeTxDbFromGRanges(mito_gff)
mito_exons <- GenomicFeatures::exonicParts(mito_txdb,
                                           linked.to.single.gene.only = TRUE)

# chloroplast
chl_gff <- subset(gff, seqnames == "ChrC")
chl_txdb <- GenomicFeatures::makeTxDbFromGRanges(chl_gff)
chl_exons <- GenomicFeatures::exonicParts(chl_txdb,
                                          linked.to.single.gene.only = TRUE)

# intergenic
# https://support.bioconductor.org/p/73648/
intergenic <- genFeatures(nuclear_txdb_with_rrna,
                          featuretype = "intergenic",
                          reduce_ranges = TRUE)$intergenic

# list of things to count
feature_list <- list(exons = nuc_exons,
                     introns = nuc_introns,
                     mitochondrion = mito_exons,
                     chloroplast = chl_exons,
                     intergenic = intergenic,
                     ribosome = ribo_gff)

feature_count_list <- lapply(feature_list,
                             ReadCountWrapper,
                             bamfile_list = bamfiles)

feature_counts <- rbindlist(feature_count_list,
                            idcol = "feature",
                            fill = TRUE,
                            use.names = TRUE)

# PLOT PER CATEGORY
feature_counts_long <- melt(feature_counts,
                            id.vars = c("feature", "gene_id", "exon_name", "feature_by", "ID", "Name", "Parent"),
                            measure.vars = names(bamfile_list),
                            variable.name = "sample_name",
                            value.name = "counts")

pd <- feature_counts_long[, .(counts = sum(counts)),
                          by = .(feature, sample_name)]
pd[, counts_per_sample := sum(counts), by = sample_name]
pd[, percentage_of_counts := counts* 100 / counts_per_sample]
pd[, sample_type := strsplit(as.character(sample_name), "_")[[1]][1],
   by = sample_name]

# plot1 <- ggplot(
#     pd,
#     aes(x = sample_name,
#         y = percentage_of_counts,
#         fill = sample_type)) +
#     scale_fill_brewer(palette = "Set1") +
#     facet_wrap(~ feature, scales = "free_y") +
#     geom_col()
# 

# INTRONIC VS EXONIC READS PER GENE 
# exonic_intronic_long <- feature_counts_long[feature %in% c("exons", "introns")]
# exonic_intronic <- dcast(exonic_intronic_long,
#                          gene_id + sample_name ~ feature,
#                          fun.aggregate = sum,
#                          value.var = "counts")
# plot2 <- ggplot(exonic_intronic, aes(x = exons, y = introns)) +
#     facet_wrap(~ sample_name) +
#     scale_x_log10() +
#     scale_y_log10() +
#     geom_smooth(method = "lm") +
#     geom_point()

# write output
fwrite(pd, fc_file)
# ggsave(filename = plot1_file,
#        plot = plot1,
#        device = cairo_pdf,
#        width = 10,
#        height = 7.5,
#        units = "in")
# ggsave(filename = plot2_file,
#        plot = plot2,
#        device = cairo_pdf,
#        width = 10,
#        height = 7.5,
#        units = "in")

# write session info
sessionInfo()
