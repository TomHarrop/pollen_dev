PlotMyGene <- function(goi) {
    my_counts <- counts(dds_lrt, normalized = TRUE)[goi,]
    my_count_dt <- data.table(sample_id = names(my_counts),
                              normalized_count = my_counts)
    stage_order <- c("UNM", "PUNM", "BCP", "TCP")
    my_count_dt[, c("stage", "plant") := tstrsplit(sample_id, "_")]
    my_count_dt[, stage := factor(stage, levels = stage_order)]
    my_count_dt[, plant := factor(plant, levels = sort(unique(plant)))]
    ggplot(my_count_dt, aes(x = stage,
                            y = normalized_count,
                            colour = plant,
                            group = 1)) +
        scale_colour_brewer(palette = "Set1") +
        scale_y_log10() +
        geom_smooth(colour = alpha("black", 0.5),
                    se = FALSE,
                    span = 0.5) +
        geom_point(size = 2,
                   alpha = 0.7,
                   shape = 16,
                   position = position_jitter(width = 0.2))
}