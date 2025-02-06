#' @title RNAmeta_peak_exon_size_plot
#' @param preliminary_analysis_data A list of data containing the peak exon size information, where each element is a list containing data for a particular analysis.
#'        The first level of the list should contain the data, and each item should include a second level with a `Sample` and `length` column.
#' @param set_group_name A vector of custom group names to replace the original sample names in the plot. The length of this vector must match the number of unique sample groups in the data.
#'        If `NULL`, the original sample names are used.
#' @returns A `ggplot` object that displays a boxplot of exon sizes (log2 scale) across samples, with facets by exon location.
#'          The boxplot colors correspond to the sample groups, and the plot includes a legend if more than one unique sample group is present.
#' @export
RNAmeta_peak_exon_size_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL){
  total_all_sites <- NULL
  for(i in 1:length(preliminary_analysis_data[[1]])){
    cct <- NULL
    all_sites <- NULL
    cct <- preliminary_analysis_data[[1]][[i]]
    cct2 = cct[[2]]
    all_sites <- data.table::as.data.table(cct2)
    all_sites[,group := IRanges::gsub("_.*", "", Sample)]
    all_sites[,Sample := IRanges::gsub(".*_", "", Sample)]
    total_all_sites <- IRanges::rbind(total_all_sites, all_sites)
  }
  if(is.null(set_group_name)){} else {
    original_group_names <- IRanges::unique(total_all_sites$Sample)
    if (length(set_group_name) != length(original_group_names)) {
      stop(paste("The number of user-defined group names does not match the number of groups in the data.",
                 "Length of set_group_name:", length(set_group_name),
                 "Length of original_group_names:", length(original_group_names)))
    }
    group_mapping <- stats::setNames(set_group_name, original_group_names)
    total_all_sites$Sample <- dplyr::recode(total_all_sites$Sample, !!!group_mapping)
  }

  if(length(IRanges::unique(total_all_sites$Sample)) > 1){
    legend_position <- "top"
  } else {
    legend_position <- "none"
  }
  p <- ggpubr::ggboxplot(total_all_sites, x = 'Sample', y = 'length', color = 'Sample',
                  facet.by = "exonLoc", yscale = "log2", width = 0.3) +
    ggplot2::scale_fill_manual(values = c("#4e62ab", "#469db4", "#87d0a8", "#ff5733", "#33ff57", "#ff33a8")) +
    ggplot2::labs(x = "", y = "log2(Length)(bp)") +
    ggplot2::theme(
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 10, colour = "#333333"),
      axis.text.y = ggplot2::element_text(vjust = 0.5),
      axis.text.x = ggplot2::element_text(size = 10, hjust = 0.5, margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(size = 10)
    ) +
    ggplot2::scale_y_continuous(trans = 'log2', breaks = scales::trans_breaks('log2', function(x) 2 ^ x)) +
    ggplot2::facet_wrap(~exonLoc, scales = 'free_x')
  return(p)
}
