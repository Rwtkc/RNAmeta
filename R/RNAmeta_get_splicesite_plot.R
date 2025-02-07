#' @title Plot Distribution of Splicing Sites
#' @param preliminary_analysis_data A list containing the data necessary for the analysis:
#'   - [2] GRanges object containing peak regions (for splicing site analysis)
#'   - [7] GRanges object with exon annotations (from TxDb or similar)
#'   - [6] Integer specifying the flank size around splicing sites
#'   - [9] Integer specifying a short flank size for the analysis
#' @param set_group_name A vector of user-defined group names for the samples. This allows the user to reassign the group names for the sample data.
#'   If this is NULL, no renaming of the groups will be done. The length of the `set_group_name` vector should match the number of unique groups in the dataset.
#' @returns A list of two `ggplot2` plots: one for the distribution at the 5' splicing site (5PSS) and one for the distribution at the 3' splicing site (3PSS).
#' @export
RNAmeta_get_splicesite_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL){
  peak_gr <- preliminary_analysis_data[[2]]
  txdb_features <- preliminary_analysis_data[[7]]
  short_flank_size <- preliminary_analysis_data[[9]]
  flank_size <- preliminary_analysis_data[[6]]
  cct = lapply(1:length(peak_gr), .get_eejunct_new, peak_gr, txdb_features$exons, short_flank_size)
  cvg_ee = data.table::rbindlist(cct)

  if(is.null(set_group_name)){} else {
    original_group_names <- IRanges::unique(cvg_ee$Sample)
    if (length(set_group_name) != length(original_group_names)) {
      stop(paste("The number of user-defined group names does not match the number of groups in the data.",
                 "Length of set_group_name:", length(set_group_name),
                 "Length of original_group_names:", length(original_group_names)))
    }
    group_mapping <- stats::setNames(set_group_name, original_group_names)
    cvg_ee$Sample <- dplyr::recode(cvg_ee$Sample, !!!group_mapping)
  }
  if(length(unique(cvg_ee$Sample)) > 1){
    legend_position <- "right"
  } else {
    legend_position <- "none"
  }

  p1 <- ggplot2::ggplot(cvg_ee[feature=="5PSS"], ggplot2::aes(x = rel_pos, colour = Sample)) +
    ggplot2::ggtitle("Distribution on 5' Splicing Sites (5PSS)") +
    ggplot2::scale_x_continuous(labels = scales::label_comma(), breaks = scales::pretty_breaks(n = 10), limits = c(-flank_size, flank_size * 2)) +
    ggplot2::scale_y_continuous(labels = scales::label_scientific(), breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::theme(axis.ticks = ggplot2::element_line(),
                   axis.line = ggplot2::element_line(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 18),
                   panel.grid.major = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 14),
                   axis.title.x = ggplot2::element_text(size = 14),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::xlab("Region around the site (0 and 5' -> 3' direction)") +
    ggplot2::ylab("Density") +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted")

  p2 <- ggplot2::ggplot(cvg_ee[feature=="3PSS"], aes(x = rel_pos, colour = Sample)) +
    ggplot2::ggtitle("Distribution on 3' Splicing Sites (3PSS)") +
    ggplot2::scale_x_continuous(labels = scales::label_comma(), breaks = scales::pretty_breaks(n = 10), limits = c(-flank_size * 2, flank_size)) +
    ggplot2::scale_y_continuous(labels = scales::label_scientific(), breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::theme(axis.ticks = ggplot2::element_line(),
                   axis.line = ggplot2::element_line(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 18),
                   panel.grid.major = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 14),
                   axis.title.x = ggplot2::element_text(size = 14),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::xlab("Region around the site (0 and 5' -> 3' direction)") +
    ggplot2::ylab("Density") +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted")

  plot_list <- list(p1, p2)
  return(plot_list)
}
