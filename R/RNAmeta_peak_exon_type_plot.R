#' @title RNAmeta_peak_exon_type_plot
#' @param preliminary_analysis_data A list of data containing the peak exon type information. The data should contain exon information, with columns such as `Sample` and `exonLoc`.
#'        Each element in the list corresponds to one set of data to be plotted.
#' @param set_group_name An optional vector of custom group names for the `Sample` column. If `NULL`, the function will use the original sample names in the data.
#'        If provided, the length of this vector must match the number of unique sample groups.
#' @returns A `ggplot` object representing a bar plot of the percentage distribution of exon locations for each sample. The bars are grouped by `Sample` and colored by `exonLoc`.
#'          The legend position will adjust based on the number of unique sample groups.
#' @export
RNAmeta_peak_exon_type_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL){
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

  all_sites <- total_all_sites
  exon_type_stat <- all_sites[, .N, by = c("exonLoc", "Sample")]
  exon_type_stat[, total := sum(N), by = c("Sample")]
  exon_type_stat[, N := N / total]

  if(is.null(set_group_name)){} else {
    original_group_names <- IRanges::unique(exon_type_stat$Sample)
    if (length(set_group_name) != length(original_group_names)) {
      stop(paste("The number of user-defined group names does not match the number of groups in the data.",
                 "Length of set_group_name:", length(set_group_name),
                 "Length of original_group_names:", length(original_group_names)))
    }
    group_mapping <- stats::setNames(set_group_name, original_group_names)
    exon_type_stat$Sample <- dplyr::recode(exon_type_stat$Sample, !!!group_mapping)
  }

  if(length(IRanges::unique(exon_type_stat$Sample)) > 1){
    legend_position <- "right"
  } else {
    legend_position <- "none"
  }

  p <- ggpubr::ggbarplot(exon_type_stat, x = 'Sample', y = 'N', fill = "exonLoc",
                  width = 0.3,
                  palette = c("#a30545", "#f26e42", "#7ecba4")) +
    ggplot2::labs(x = "", y = "Percentage") +
    ggplot2::theme(legend.position = legend_position,
          legend.title = ggplot2::element_text(size = 10),
          text = ggplot2::element_text(color = "#333333", size = 8),
          legend.text = ggplot2::element_text(size = 12),
          axis.title.y = ggplot2::element_text(size = 12),
          axis.text.x = ggplot2::element_text(size = 12),
          axis.text.y = ggplot2::element_text(size = 10),
          axis.ticks = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line(color = "gray90"),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(color = "#333333")) + ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 100))
  return(p)
}
