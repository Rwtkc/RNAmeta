#' @title RNA Meta Peak Gene Size Plot
#' @param preliminary_analysis_data A list containing the preliminary analysis data. It should include the transcript length (`tx_len`) information extracted from the dataset.
#' @param set_group_name A vector of user-defined group names to rename the sample groups in the dataset.
#'   The length of `set_group_name` must match the number of unique `Sample` values in the data.
#' @returns A `ggplot` object representing a box plot of transcript lengths (`tx_len`) across different sample groups.
#'   - The x-axis represents different sample groups.
#'   - The y-axis represents transcript length in log2 scale (`log2(txLength bp)`).
#'   - The box plot shows the distribution of transcript lengths within each sample group.
#'   - The legend placement is determined based on the number of unique samples.
#'
#' @export
#'
RNAmeta_peak_gene_size_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL){
  total_all_sites <- NULL
  for(i in 1:length(preliminary_analysis_data[[1]])){
    all_sites <- NULL
    cct <- NULL
    cct <- preliminary_analysis_data[[1]][[i]]
    cct7 = cct[5]
    all_sites = data.table::rbindlist(cct7)
    all_sites[,group := gsub("_.*", "", Sample)]
    all_sites[,Sample := gsub(".*_", "", Sample)]
    data.table::as.data.table(all_sites)
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
  p <- ggpubr::ggboxplot(total_all_sites, x = 'Sample', y = 'tx_len', color = "Sample", yscale ="log2", width = 0.3) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(),
      axis.ticks.y = ggplot2::element_line(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 10, colour = "#333333"),
      axis.text.y = ggplot2::element_text(vjust = 0.5),
      axis.text.x = ggplot2::element_text(size = 10,hjust = 0.5,margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(r = 10))) +
    ggplot2::labs(y = "log2(txLenght bp)")
  return(p)
}
