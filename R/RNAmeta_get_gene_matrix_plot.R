#' @title Generate a Venn Diagram of Gene Matrix
#' @param preliminary_analysis_data A list containing analysis data. The first element should be a list of data frames or similar structures where each element has a 'gene_name' column.
#' @param set_group_name A character vector of group names. If provided, it will be used to label the Venn diagram categories. The length of this vector should match the number of elements in `preliminary_analysis_data[[1]]`.
#'   If `NULL`, default group names ("Group1", "Group2", "Group3", "Group4", "Group5") will be used.
#' @returns A Venn diagram object (from the `VennDiagram` package) representing the gene overlap between the BED files.
#' @export
RNAmeta_get_gene_matrix_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL){
  cct <- preliminary_analysis_data[[1]]
  cct_file_length <- length(cct)
  if (cct_file_length < 2 || cct_file_length > 5) {
    stop("Error: The number of BED files required for the gene matrix should be between 2 and 5. Current number of files: ", cct_file_length)
  }
  all_genes_data_list <- list()
  for(i in 1:cct_file_length){
    all_genes <- NULL
    all_genes <- as.data.frame(unique(na.omit(cct[[i]][[3]]$gene_name)), stringsAsFactors = FALSE)
    all_genes_data_list[[i]] <- all_genes
  }
  venn_color <- c("#FF0000","#ff8080","#a1c530","#6676a0","#33b39f")
  venn_group_name <- c("Group1", "Group2","Group3","Group4","Group5")
  if(is.null(set_group_name)){
    set_group_name <- venn_group_name
  } else {
    if (length(set_group_name) > length(cct)) {
      set_group_name <- set_group_name[1:length(cct)]
    } else if (length(set_group_name) < length(cct)) {
      set_group_name <- venn_group_name
    }
  }
  gene_venn_data <- lapply(all_genes_data_list, function(df) df[, 1])
  veen_plot <- VennDiagram::venn.diagram(
    x = gene_venn_data,
    disable.logging = T,
    category.names = set_group_name[1:length(all_genes_data_list)],
    fill = venn_color[1:length(all_genes_data_list)],
    filename = NULL,
    cat.col = venn_color[1:length(all_genes_data_list)],
    col = "white",
    cat.cex = 0.8,
    margin = 0.1
  )
  return(veen_plot)
}
