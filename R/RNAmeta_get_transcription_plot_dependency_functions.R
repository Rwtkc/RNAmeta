.get_feature_boundary_coverage_new <- function(index = NULL, query_regions = NULL, feature_coords = NULL, flank_size = NULL, boundary_type = NULL,  sample_n = NULL) {
  if(length(feature_coords) == 0){
    return (data.table::data.table())
  }
  feature_coords = unlist(feature_coords)
  if (sample_n > 0 && sample_n < length(feature_coords)) {
    feature_coords <- sort(feature_coords[sample(length(feature_coords), sample_n)])
  }

  if (boundary_type == 'Transcription_start_site') {
    flanks <- GenomicRanges::flank(x = feature_coords,
                                   width = flank_size,
                                   start = TRUE,
                                   both = TRUE)
  } else if (boundary_type == 'Transcription_end_site') {
    flanks <- GenomicRanges::flank(x = feature_coords,
                                   width = flank_size,
                                   start = FALSE,
                                   both = TRUE)
  } else {
    stop ("please indicate either Transcription_start_site or Transcription_end_site for boundary type\n")
  }

  sample_name <- names(query_regions[index])
  peak_table <- data.table::as.data.table(query_regions[[index]])
  peak_table[, index:=.I]
  hits <- GenomicFeatures::mapToTranscripts(query_regions[[index]], flanks, ignore.strand = FALSE)
  hits <- data.table::as.data.table(hits)
  setnames(hits, c('seqnames', 'xHits'), c('transcriptID', 'index'))
  peak_table <- peak_table[hits, on ='index']
  peak_table[, c("index", "transcriptsHits") := NULL]
  peak_table[, Sample := sample_name]
  peak_table[, feature := boundary_type]
  peak_table[, rel_pos := i.start - flank_size]
  peak_table[, flank_size := flank_size]
  peak_table[, xmin := -flank_size]
  peak_table[, xmax := flank_size - 1]
  return(peak_table)
}

.score_matrix_RNA <-function(target = NULL, windows = NULL, strand_aware = NULL, weight_col = NULL){
  if(is.null(weight_col)){
    target.rle = SummarizedExperiment::coverage(target)
  }else{
    if(!weight_col %in% names(GenomicRanges::mcols(target)) ){
      stop("provided column 'weight_col' does not exist in target\n")
    }
    target.rle = SummarizedExperiment::coverage(target, weight = weight_col)
  }
  if( length(unique(rtracklayer::width(windows))) > 1 ){
    stop("width of 'windows' are not equal, provide 'windows' with equal widths")
  }
  if( any(rtracklayer::width(windows) == 1) ){
    stop("provide 'windows' with widths greater than 1")
  }
  windows_length = length(windows)
  windows = genomation::constrainRanges(target.rle, windows)
  if(length(windows) != windows_length){
    warning(paste0(windows_length - length(windows),
                   " windows fall off the target"))
  }

  chrs = sort(dplyr::intersect(names(target.rle), as.character(unique(SummarizedExperiment::seqnames(windows)))))
  my_views = IRanges::Views(target.rle[chrs], as(windows,"IntegerRangesList")[chrs])
  mat = lapply(my_views, function(x) t(IRanges::viewApply(x, as.vector)) )
  mat = do.call("rbind", mat)

  r_list = split(GenomicRanges::mcols(windows)[, "X_rank"], as.vector(SummarizedExperiment::seqnames(windows)))
  r_list = r_list[order(names(r_list))]
  ranks = do.call("c", r_list)
  rownames(mat) = ranks

  if(strand_aware == TRUE){
    orig_rows = windows[BiocGenerics::strand(windows) == '-', ]$X_rank
    mat[rownames(mat) %in% orig_rows, ] = mat[rownames(mat) %in%
                                                orig_rows, ncol(mat):1]
  }
  mat = mat[order(ranks), ]
  return(mat)
}

.peak_heatmap <- function(tag_matrix_5 = NULL, tag_matrix_3 = NULL, color = "red") {
  ii <- order(rowSums(tag_matrix_5))
  tag_matrix_5 <- tag_matrix_5[ii, ]
  if(length(which(rowSums(tag_matrix_5) > 0)) <= 1) {
    rowN_5 <- ifelse(nrow(tag_matrix_5) >= 1000, 1000, nrow(tag_matrix_5))
    tag_matrix_5 <- tag_matrix_5[1:rowN_5, ]
    tag_matrix_5[1, 1] = 0.5
  } else {
    tag_matrix_5 <- tag_matrix_5[rowSums(tag_matrix_5) > 0, ]
  }

  xlim_5 <- seq(-ncol(tag_matrix_5) / 2 + 1, ncol(tag_matrix_5) / 2)


  ii <- order(rowSums(tag_matrix_3))
  tag_matrix_3 <- tag_matrix_3[ii, ]
  if(length(which(rowSums(tag_matrix_3) > 0)) <= 1){
    rowN_3 <- ifelse(nrow(tag_matrix_3) >= 1000, 1000, nrow(tag_matrix_3))
    tag_matrix_3 <- tag_matrix_3[1:rowN_3, ]
    tag_matrix_3[1, 1] = 0.5
  } else {
    tag_matrix_3 <- tag_matrix_3[rowSums(tag_matrix_3) > 0, ]
  }
  xlim_3 <- seq(-ncol(tag_matrix_3) / 2 + 1,ncol(tag_matrix_3) / 2)
  cols <- grDevices::colorRampPalette(c("white", color))(200)
  peak_heatmap_data <- list(xlim_3, xlim_5, tag_matrix_3, tag_matrix_5)
  return(peak_heatmap_data)
}

.get_feature_boundary_heatmap <- function(index = NULL, query_regions = NULL, txdb_features = NULL, gene_max_tx_lengths = NULL, flank_size = NULL) {
  sn = names(query_regions[index])
  qregion = query_regions[[index]]
  qregion <- GenomicRanges::resize(qregion, 15, fix="center")
  flanks_5 <- GenomicRanges::flank(x = txdb_features$five_utrs,
                                   width = flank_size,
                                   start = TRUE,
                                   both = TRUE)
  flanks_3 <- GenomicRanges::flank(x = txdb_features$three_utrs,
                                   width = flank_size,
                                   start = FALSE,
                                   both = TRUE)
  flanks_5 <- flanks_5[names(flanks_5) %in% gene_max_tx_lengths$tx_name]
  flanks_3 <- flanks_3[names(flanks_3) %in% gene_max_tx_lengths$tx_name]
  flanks_5 <- BiocGenerics::unlist(flanks_5)
  flanks_3 <- BiocGenerics::unlist(flanks_3)
  sm_5 <- suppressWarnings(.score_matrix_RNA(target = qregion,
                                             windows = flanks_5,
                                             strand_aware = TRUE))
  sm_3 <- suppressWarnings(.score_matrix_RNA(target = qregion,
                                             windows = flanks_3,
                                             strand_aware = TRUE))

  return(.peak_heatmap(tag_matrix_5 = sm_5, tag_matrix_3 = sm_3))
}

.get_transcription_data <- function(preliminary_analysis_data = NULL){
  txdb_features <- preliminary_analysis_data[[7]]
  peak_gr <- preliminary_analysis_data[[2]]
  flank_size <- preliminary_analysis_data[[6]]
  gene_max_tx_lengths <- preliminary_analysis_data[[8]]
  cct3 = lapply(1:length(peak_gr), .get_feature_boundary_coverage_new, query_regions = peak_gr, feature_coords = txdb_features$cds, flank_size = flank_size, boundary_type = 'Transcription_start_site', sample_n = 10000)
  cct4 = lapply(1:length(peak_gr), .get_feature_boundary_coverage_new , query_regions = peak_gr, feature_coords = txdb_features$cds, flank_size = flank_size, boundary_type = 'Transcription_end_site', sample_n = 10000)
  cvg_list = data.table::rbindlist(c(cct3,cct4))
  peak_heatmap_data = lapply(1:length(peak_gr), .get_feature_boundary_heatmap, peak_gr, txdb_features, gene_max_tx_lengths, flank_size)
  cvg_list_and_peak_heatmap_data <- list(cvg_list, peak_heatmap_data[[1]])
  return(cvg_list_and_peak_heatmap_data)
}
