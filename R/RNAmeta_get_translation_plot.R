#' @title Get Translation Data
#' @param preliminary_analysis_data A list containing multiple elements:
#'   - The 7th element (`txdb_features`) contains genomic feature data (e.g., exons).
#'   - The 2nd element (`peak_gr`) contains the genomic regions of interest (e.g., peaks).
#'   - The 9th element (`short_flank_size`) represents the flank size used for analysis.
#'   - The 10th element (`tx_length_file`) contains transcript length information for the genes.
#' @returns A list with two elements:
#'   1. A data table containing the coverage information for both translation start sites and end sites (`tsl_list`).
#'   2. A list of peak heatmap data generated from the translation-related regions (`peak_heatmap_data`).
#' @export
get_translation_data <- function(preliminary_analysis_data = NULL){
  txdb_features <- preliminary_analysis_data[[7]]
  peak_gr <- preliminary_analysis_data[[2]]
  short_flank_size <- preliminary_analysis_data[[9]]
  tx_length_file <- preliminary_analysis_data[[10]]
  cct_tran = lapply(1:length(peak_gr), .get_translation_boundary_coverage_new, peak_gr, tx_length_file, txdb_features$exons, short_flank_size)
  tsl_list <- data.table::rbindlist(cct_tran)
  cct_tran = lapply(1:length(peak_gr), .get_translation_boundary_heatmap, peak_gr, tx_length_file, txdb_features$exons, short_flank_size)
  tsl_list_and_cct_tran <- list(tsl_list, cct_tran[[1]])
  return(tsl_list_and_cct_tran)
}

#' @title Generate Translation Plot
#' @param translation_data A list containing:
#'   - `tsl_list`: A data table with translation start and end site coverage.
#'   - `peak_heatmap_data`: A list containing translation heatmap data.
#' @param set_group_name An optional vector of custom group names to replace the default sample names in the data.
#'        If `NULL`, the function will use the original sample names.
#' @returns A list of four ggplot objects:
#'   1. Translation Start Site (TSS) Density Plot.
#'   2. Translation End Site (TES) Density Plot.
#'   3. Heatmap for the TSS region.
#'   4. Heatmap for the TES region.
#' @export
RNAmeta_get_translation_plot <- function(translation_data = NULL, set_group_name = NULL){
  tsl_list <- translation_data[[1]]
  peak_heatmap_data <- translation_data[[2]]
  cct3 <- tsl_list[feature=="Translation_start_site"]
  xmin <- unique(cct3$xmin)
  xmax <- unique(cct3$xmax)
  flank_size_3 <- unique(cct3$flank_size)
  if(is.null(set_group_name)){} else {
    original_group_names <- IRanges::unique(cct3$Sample)
    if (length(set_group_name) != length(original_group_names)) {
      stop(paste("The number of user-defined group names does not match the number of groups in the data.",
                 "Length of set_group_name:", length(set_group_name),
                 "Length of original_group_names:", length(original_group_names)))
    }
    group_mapping <- stats::setNames(set_group_name, original_group_names)
    cct3$Sample <- dplyr::recode(cct3$Sample, !!!group_mapping)
  }
  if(length(unique(cct3$Sample)) > 1){
    legend_position <- "right"
  } else {
    legend_position <- "none"
  }
  p1 <- ggplot2::ggplot(cct3, ggplot2::aes(x = rel_pos, colour = Sample)) +
    ggplot2::ggtitle("Distribution on Translation start site") +
    ggplot2::scale_x_continuous(labels = scales::label_comma(), breaks = scales::pretty_breaks(n = 10), limits = c(xmin, xmax)) +
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
  #----------------------------------------------------
  cct4 <- tsl_list[feature=="Translation_end_site"]
  xmin2 <- unique(cct4$xmin)
  xmax2 <- unique(cct4$xmax)
  flank_size_4 <- unique(cct4$flank_size)
  if(is.null(set_group_name)){} else {
    original_group_names <- IRanges::unique(cct4$Sample)
    if (length(set_group_name) != length(original_group_names)) {
      stop(paste("The number of user-defined group names does not match the number of groups in the data.",
                 "Length of set_group_name:", length(set_group_name),
                 "Length of original_group_names:", length(original_group_names)))
    }
    group_mapping <- stats::setNames(set_group_name, original_group_names)
    cct4$Sample <- dplyr::recode(cct4$Sample, !!!group_mapping)
  }
  if(length(unique(cct4$Sample)) > 1){
    legend_position <- "right"
  } else {
    legend_position <- "none"
  }
  p2 <- ggplot2::ggplot(cct4, ggplot2::aes(x = rel_pos, colour = Sample)) +
    ggplot2::ggtitle("Distribution on Translation end site") +
    ggplot2::scale_x_continuous(labels = scales::label_comma(), breaks = scales::pretty_breaks(n = 10), limits = c(xmin2, xmax2)) +
    ggplot2::scale_y_continuous(labels = scales::label_scientific(), breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::theme(axis.ticks = ggplot2::element_line(),
                   axis.line = ggplot2::element_line(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 18),
                   axis.title.y = ggplot2::element_text(size = 14),
                   axis.title.x = ggplot2::element_text(size = 14),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::xlab("Region around the site (0 and 5' -> 3' direction)") +
    ggplot2::ylab("Density") +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
  #----------------------------------------------------
  xlim_5 <- peak_heatmap_data[[2]]
  xlim_3 <- peak_heatmap_data[[1]]
  tag_matrix_5 <- peak_heatmap_data[[4]]
  tag_matrix_3 <- peak_heatmap_data[[3]]
  min_xlim_5 <- min(xlim_5)
  max_xlim_5 <- max(xlim_5)
  rows_5 <- nrow(tag_matrix_5)
  cols_5 <- ncol(tag_matrix_5)

  indices_5 <- S4Vectors::expand.grid(row = 1:rows_5, col = 1:cols_5)
  values_5 <- apply(indices_5, 1, function(idx) tag_matrix_5[idx["row"], idx["col"]])
  result_data_5 <- data.frame(row_number = indices_5$row,
                              col_number = indices_5$col,
                              value = values_5)
  data_column_5 <- integer(0)
  for (i in min_xlim_5 : max_xlim_5 ) {
    data_column_5 <- c(data_column_5, rep(i, rows_5))
  }
  result_data_5$col_number <- data_column_5
  p3 <- ggplot2::ggplot(result_data_5, ggplot2::aes(x = col_number, y = row_number, fill = value)) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::scale_fill_gradientn(colors = c("#420046", "#fedc00")) +
    ggplot2::scale_x_continuous(expand = c(0.003, 0), breaks = pretty(result_data_5$col_number, n = 8)) +
    ggplot2::scale_y_continuous(expand = c(0.008, 0), breaks = pretty(result_data_5$row_number, n = 10)) +
    ggplot2::labs(title = "Transcription.Tss.Heatmap") +
    ggplot2::theme_minimal() + ggplot2::theme(
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(5.5, 14, 5.5, 14), "points"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  #----------------------------------------------------
  min_xlim_3 <- min(xlim_3)
  max_xlim_3 <- max(xlim_3)
  rows_3 <- nrow(tag_matrix_3)
  cols_3 <- ncol(tag_matrix_3)
  indices_3 <- S4Vectors::expand.grid(row = 1:rows_3, col = 1:cols_3)
  values_3 <- apply(indices_3, 1, function(idx) tag_matrix_3[idx["row"], idx["col"]])
  result_data_3<- data.frame(row_number = indices_3$row,
                             col_number = indices_3$col,
                             value = values_3)
  data_column_3 <- integer(0)
  for (i in min_xlim_3 : max_xlim_3 ) {
    data_column_3 <- c(data_column_3, rep(i, rows_3))
  }
  result_data_3$col_number <- data_column_3

  p4 <- ggplot2::ggplot(result_data_3, aes(x = col_number, y = row_number, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = c("#420046", "#fedc00")) +
    ggplot2::scale_x_continuous(expand = c(0.003, 0), breaks = pretty(result_data_3$col_number, n = 8)) +
    ggplot2::scale_y_continuous(expand = c(0.008, 0), breaks = pretty(result_data_3$row_number, n = 10)) +
    ggplot2::labs(title = "Translation.Tes.Heatmap") +
    ggplot2::theme_minimal() + ggplot2::theme(
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(5.5, 14, 5.5, 14), "points"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  #----------------------------------------------------
  plot_list <- list(p1, p2, p3, p4)
  return(plot_list)
}


