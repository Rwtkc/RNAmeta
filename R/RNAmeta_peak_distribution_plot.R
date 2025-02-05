#' @title Plot Peak Distribution
#'
#' @param preliminary_analysis_data A list containing the preliminary analysis data, typically including genomic features and their frequencies across different samples.
#' @param set_group_name A vector of group names to be used for color grouping in the plot. The number of group names should match the number of unique samples in `all_sites`.
#' @param features_to_remove A vector of feature names (e.g., "Promoter", "UTR5", etc.) that should be excluded from the analysis and visualization.
#'
#' @returns A `ggplot` object representing the bar plot of peak distribution across features. The plot shows the frequency of each feature for different samples, with a manual color scheme applied to each group.
#' @export
#'
RNAmeta_peak_distribution_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL, features_to_remove = NULL){
  all_sites_list <- list()
  my_colors <- c("#de7a61", "#3f3f5b", "#81b29d", "#ffb5b0", "#d3a6c6", "#f5d300")
  for(i in 1:length(preliminary_analysis_data[[1]])){
    cct <- preliminary_analysis_data[[1]][[i]]
    cct1 = cct[[1]]
    all_sites <- list(cct1)
    desired_order <- c("Promoter", "UTR5", "Start Codon", "CDS", "Stop Codon", "UTR3", "Intron", "Intergenic")
    all_sites[[1]]$Feature <- factor(all_sites[[1]]$Feature, levels = desired_order)
    all_sites[[1]] <- all_sites[[1]][order(all_sites[[1]]$Feature), ]
    all_sites_list <- IRanges::append(all_sites_list, all_sites)
  }
  all_sites <- NULL
  for(i in 1:length(all_sites_list)){
    all_sites <- IRanges::rbind(all_sites,all_sites_list[[i]])
  }
  all_sites <- data.table::as.data.table(all_sites)
  features <- c("Promoter", "UTR5", "Start Codon", "CDS", "Stop Codon", "UTR3", "Intron", "Intergenic")
  features <- setdiff(features, features_to_remove)
  all_sites <- AnnotationDbi::subset(all_sites, Feature %in% features)
  if (length(IRanges::unique(all_sites$Sample)) == 1) {
    legend_position <- "none"
  } else {
    legend_position <- "right"
  }
  if(is.null(set_group_name)){
    set_group_name <- unique(all_sites$Sample)
  } else if(length(set_group_name) != length(unique(all_sites$Sample))){
    stop(paste("The number of group names does not match the number of unique samples.",
               "Length of set_group_name:", length(set_group_name),
               "Length of unique(all_sites$Sample):", length(unique(all_sites$Sample))))
  }
  p <- ggpubr::ggbarplot(all_sites, x = "Feature", y = "Frequency",
                         color = "black",
                         fill = "Sample",
                         position = ggplot2::position_dodge(0.8),
                         label = TRUE,
                         lab.pos = "out",
                         lab.vjust = -0.5) +
    ggplot2::labs(y = "Frequency", x = "") +
    ggplot2::theme(axis.ticks.y = ggplot2::element_line(color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_line(color = "black"),
                   axis.line.x = ggplot2::element_line(color = "black"),
                   legend.position = legend_position,
                   axis.text.x = ggplot2::element_text(angle = 0, vjust = -1),
                   axis.text.y = ggplot2::element_text(angle = 0, hjust = 0.1),
                   axis.title.y = ggplot2::element_text(vjust = 2.5, size = 15),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(labels = function(y) paste0(round(y / 1000, 1), "k"), expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::scale_fill_manual(values = my_colors, labels = set_group_name)
  return(p)
}
