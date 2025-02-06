#' @title RNAmeta_peak_exon_number_plot
#' @param preliminary_analysis_data A list containing exon-related data. The first level of the list contains data for different analyses, and each analysis includes a sublist. The specific data needed is located at the 5th index in each sublist, containing the `Sample` and `nexon` (number of exons) information.
#' @param set_group_name An optional vector of custom group names to replace the original sample names in the data. If `NULL`, the original sample names are retained. The length of this vector must match the number of unique sample groups.
#' @returns A `ggplot` object representing a boxplot showing the distribution of the number of exons across different sample groups. The y-axis represents the number of exons, and the boxplots are grouped by `Sample`.
#' @export
RNAmeta_peak_exon_number_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL){
  total_all_sites <- NULL
  for(i in 1:length(preliminary_analysis_data[[1]])){
    all_sites <- NULL
    cct <- NULL
    cct <- preliminary_analysis_data[[1]][[i]]
    cct7 = cct[5]
    all_sites = data.table::rbindlist(cct7)
    all_sites[, group := gsub("_.*", "", Sample)]
    all_sites[, Sample := gsub(".*_", "", Sample)]
    data.table::as.data.table(all_sites)
    total_all_sites <- IRanges::rbind(total_all_sites, all_sites)
  }
  all_sites <- total_all_sites
  if(is.null(set_group_name)){} else {
    original_group_names <- IRanges::unique(all_sites$Sample)
    if (length(set_group_name) != length(original_group_names)) {
      stop(paste("The number of user-defined group names does not match the number of groups in the data.",
                 "Length of set_group_name:", length(set_group_name),
                 "Length of original_group_names:", length(original_group_names)))
    }
    group_mapping <- stats::setNames(set_group_name, original_group_names)
    all_sites$Sample <- dplyr::recode(all_sites$Sample, !!!group_mapping)
  }

  if(length(IRanges::unique(all_sites$Sample)) > 1){
    legend_position <- "top"
  } else {
    legend_position <- "none"
  }
  p <- ggpubr::ggboxplot(all_sites, x = 'Sample', y = 'nexon', width = 0.3, color = "Sample") +
    ggplot2::labs(x = "", y = "Number of exons") +
    ggplot2::scale_fill_manual(values = c("#c21d2e", "#7ecba4", "#ff6347", "#4682b4", "#32cd32", "#ffd700")) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(),
      axis.ticks.y = ggplot2::element_line(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 10, colour = "#333333"),
      axis.text.y = ggplot2::element_text(vjust = 0.5),
      axis.text.x = ggplot2::element_text(size = 10, hjust = 0.5, margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(r = 10))
    )
  return(p)
}
