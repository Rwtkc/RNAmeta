#' @title RNAmeta_get_transcription_plot
#' @param transcription_data A list containing two elements:
#'   - The first element is a data frame or data table containing information on transcription coverage (`cvg_list`).
#'   - The second element is a list with data related to peak heatmap data for the transcriptional regions (i.e., `peak_heatmap_data`).
#' @param set_group_name An optional vector of custom group names to replace the default sample names in the data. The length of this vector must match the number of unique sample groups in the data.
#'        If `NULL`, the function will use the original sample names.
#' @returns A list of ggplot objects. The list contains the following plots:
#'   1. A density plot of transcription coverage around the transcription start site (TSS).
#'   2. A density plot of transcription coverage around the transcription end site (TES).
#'   3. A heatmap of transcription signal for the TSS region.
#'   4. A heatmap of transcription signal for the TES region.
#' @export
RNAmeta_get_transcription_plot <- function(transcription_data = NULL, set_group_name = NULL){
  cvg_list <- transcription_data[[1]]
  peak_heatmap_data <- transcription_data[[2]]

  col <- grDevices::colorRampPalette(c("#420046", "#fedc00"))(200)
  if(nrow(cvg_list) > 0){
    cct3 <- cvg_list[feature == "Transcription_start_site"]
    cct4 <- cvg_list[feature == "Transcription_end_site"]
    xmin <- IRanges::unique(cct3$xmin)
    xmax <- IRanges::unique(cct3$xmax)
    flank_size_3 <- IRanges::unique(cct3$flankSize)

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

    p1 <- suppressWarnings(ggplot2::ggplot(cct3, ggplot2::aes(x = rel_pos, colour = Sample)) +
                             ggplot2::ggtitle("Distribution on transcription start site") +
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
                             ggplot2::geom_vline(xintercept = 0, linetype = "dotted"))
    #----------------------------------------------------
    xmin2 <- IRanges::unique(cct4$xmin)
    xmax2 <- IRanges::unique(cct4$xmax)
    flank_size_4 <- IRanges::unique(cct4$flankSize)
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
    if(length(IRanges::unique(cct4$Sample)) > 1){
      legend_position <- "right"
    } else {
      legend_position <- "none"
    }
    p2 <- suppressWarnings(ggplot2::ggplot(cct4, ggplot2::aes(x = rel_pos, colour = Sample)) +
                             ggplot2::ggtitle("Distribution on transcription end site") +
                             ggplot2::scale_x_continuous(labels = scales::label_comma(), breaks = scales::pretty_breaks(n = 10), limits = c(xmin2, xmax2)) +
                             ggplot2::scale_y_continuous(labels = scales::label_scientific(), breaks = scales::pretty_breaks(n = 8)) +
                             ggplot2::theme(axis.ticks = ggplot2::element_line(),
                                            axis.line = ggplot2::element_line(),
                                            plot.title = ggplot2::element_text(hjust = 0.5,size = 18),
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
                             ggplot2::geom_vline(xintercept = 0, linetype = "dotted"))
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
    p3 <- suppressWarnings(ggplot2::ggplot(result_data_5, ggplot2::aes(x = col_number, y = row_number, fill = value)) +
                             ggplot2::geom_tile() +
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
                               plot.title = ggplot2::element_text(hjust = 0.5))
    )
    #----------------------------------------------------
    min_xlim_3 <- min(xlim_3)
    max_xlim_3 <- max(xlim_3)
    rows_3 <- nrow(tag_matrix_3)
    cols_3 <- ncol(tag_matrix_3)
    indices_3 <- S4Vectors::expand.grid(row = 1:rows_3, col = 1:cols_3)
    values_3 <- apply(indices_3, 1, function(idx) tag_matrix_3[idx["row"], idx["col"]])
    result_data_3 <- data.frame(row_number = indices_3$row,
                                col_number = indices_3$col,
                                value = values_3)
    data_column_3 <- integer(0)
    for (i in min_xlim_3 : max_xlim_3 ) {
      data_column_3 <- c(data_column_3, rep(i, rows_3))
    }
    result_data_3$col_number <- data_column_3
    p4 <- suppressWarnings(ggplot2::ggplot(result_data_3, ggplot2::aes(x = col_number, y = row_number, fill = value)) +
                             ggplot2::geom_tile() +
                             ggplot2::scale_fill_gradientn(colors = c("#420046", "#fedc00")) +
                             ggplot2::scale_x_continuous(expand = c(0.003, 0), breaks = pretty(result_data_3$col_number, n = 8)) +
                             ggplot2::scale_y_continuous(expand = c(0.008, 0), breaks = pretty(result_data_3$row_number, n = 10)) +
                             ggplot2::labs(title = "Transcription.Tes.Heatmap") +
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
                             ))
    plot_list <- list(p1, p2, p3, p4)
    return(plot_list)
  }
}
