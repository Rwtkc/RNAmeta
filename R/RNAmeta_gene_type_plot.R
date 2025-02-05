#' @title Gene Type Plot
#'
#' @param preliminary_analysis_data A list containing the preliminary analysis data. This typically includes genomic ranges, GFF file, and other results.
#'   - The second element should be genomic ranges (`peak_gr`).
#'   - The third element should be a GFF file (`gff_file`).
#'
#' @param set_group_name A vector of user-defined group names to map to the unique `Sample` identifiers in the `all_sites` data.
#'   The length of `set_group_name` must match the number of unique `Sample` values in the data.
#'
#' @returns A `ggplot` object representing a bar plot of gene types in the dataset. The plot shows the count of each gene type across samples,
#'   with bars colored based on the sample groups. The labels indicate the count values, and the position of the legend depends on the number
#'   of unique samples (if more than one group is present, the legend is shown).
#'
#' @export
#'
RNAmeta_gene_type_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL){
  peak_gr <- preliminary_analysis_data[[2]]
  gff_file <- preliminary_analysis_data[[3]]
  cct_res <- list()

  for(i in 1:length(preliminary_analysis_data[[1]])){
    cct <- NULL
    cct = IRanges::lapply(1:length(peak_gr[i]), .batch_peak_type, peak_gr[i], gff_file)
    cct_res <- IRanges::append(cct, cct_res)
  }
  all_sites <- NULL
  for(i in 1:length(cct_res)){
    all_sites <- IRanges::rbind(cct_res[[i]], all_sites)
  }
  all_sites <- data.table::as.data.table(all_sites)

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
    legendPosition <- "right"
  } else {
    legendPosition <- "none"
  }

  p <- suppressWarnings(ggpubr::ggbarplot(all_sites, x = "V1", y = "N", fill = "Sample",
                                          label = TRUE,
                                          lab.pos = "out",
                                          lab.hjust = -0.7,
                                          lab.vjust = 0.5,
                                          lab.size = 4,
                                          sort.by.groups = TRUE,
                                          rotate = TRUE,
                                          position = ggplot2::position_dodge(0.8)) +
                          ggplot2::labs(y = "", x = "") +
                          ggplot2::theme(legend.position = legendPosition,
                                         axis.ticks.y = ggplot2::element_blank(),
                                         axis.ticks.x = ggplot2::element_blank(),
                                         panel.grid.major.x = ggplot2::element_line(color = "grey", size = 0.45),
                                         axis.line.x = ggplot2::element_blank(),
                                         axis.text.x = ggplot2::element_text(color = "#333333", size = ggplot2::rel(1), angle=0, margin = ggplot2::margin(t = 10, b = 0)),
                                         axis.text.y = ggplot2::element_text(color = "#333333", size = ggplot2::rel(1), margin = ggplot2::margin(r = 5))
                          ) +
                          ggplot2::scale_y_continuous(
                            labels = function(y) paste0(round(y / 1000, 1), "k"),
                            expand = ggplot2::expansion(mult = c(0, 0.2))
                          ))
  return(p)
}
