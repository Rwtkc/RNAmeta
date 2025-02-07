.get_translation_boundary_coverage_new <- function(index, query_regions, tx_length_file, feature_coords, flank_size = NULL) {
  tx_length_for_translation_site <- data.table::copy(tx_length_file)
  tx_length_for_translation_site <- tx_length_for_translation_site[utr5_len > flank_size & utr3_len > flank_size & cds_len > flank_size, .(tx_name, utr5_len, cds_len, utr3_len)]
  tss_gene_ranges <- tx_length_for_translation_site[, .(seqnames = tx_name, start = utr5_len + 1 - flank_size, end = utr5_len + 1 + flank_size, strand = "*", tx_name = tx_name)]
  tes_gene_ranges <- tx_length_for_translation_site[, .(seqnames = tx_name, start = utr5_len + cds_len - flank_size, end = utr5_len + cds_len + flank_size, strand = "*", tx_name = tx_name)]
  tss_gene_ranges <- GenomicRanges::makeGRangesFromDataFrame(tss_gene_ranges)
  tes_gene_ranges <- GenomicRanges::makeGRangesFromDataFrame(tes_gene_ranges)
  tss_gene_ranges <- GenomicRanges::split(tss_gene_ranges, GenomicRanges::seqnames(tss_gene_ranges))
  tes_gene_ranges <- GenomicRanges::split(tes_gene_ranges, GenomicRanges::seqnames(tes_gene_ranges))
  qregion <- query_regions[[index]]
  qregion <- qregion[as.character(GenomicRanges::seqnames(qregion)) %in% as.character(GenomeInfoDb::seqlevels(feature_coords))]
  GenomeInfoDb::seqlevels(qregion) <- intersect(GenomeInfoDb::seqlevels(qregion), GenomeInfoDb::seqlevels(feature_coords))
  peak_gene_ranges <- GenomicFeatures::mapToTranscripts(qregion, feature_coords)
  results_list <- list()
  if(length(peak_gene_ranges) == 0){
    results_list$tss <- data.table::data.table()
    results_list$tes <- data.table::data.table()
    results_list <- data.table::rbindlist(results_list)
    return(results_list)
  }
  sample_name <- names(query_regions[index])
  hits_tss <- GenomicFeatures::mapToTranscripts(peak_gene_ranges, tss_gene_ranges, ignore.strand = FALSE)
  peak_table <- data.table::as.data.table(query_regions[[index]])
  peak_table[, index := .I]
  hits_tss <- data.table::as.data.table(hits_tss)
  data.table::setnames(hits_tss, c('seqnames','xHits'), c('transcriptID','index'))
  peak_table_tss <- peak_table[hits_tss, on = 'index']
  peak_table_tss[, c("index","transcriptsHits") := NULL]
  peak_table_tss[, Sample := sample_name]
  peak_table_tss[, feature := 'Translation_start_site']
  peak_table_tss[, rel_pos := i.start - flank_size]
  peak_table_tss[, flank_size := flank_size]
  hits_tes <- GenomicFeatures::mapToTranscripts(peak_gene_ranges, tes_gene_ranges, ignore.strand = FALSE)
  hits_tes <- data.table::as.data.table(hits_tes)
  data.table::setnames(hits_tes, c('seqnames','xHits'), c('transcriptID','index'))
  peak_table_tes <- peak_table[hits_tes, on = 'index']
  peak_table_tes[, c("index","transcriptsHits") := NULL]
  peak_table_tes[, Sample := sample_name]
  peak_table_tes[, feature := 'Translation_end_site']
  peak_table_tes[, rel_pos := i.start - flank_size]
  peak_table_tes[, flank_size := flank_size]
  peak_table_tss[, xmin := -flank_size]
  peak_table_tss[, xmax := flank_size]
  peak_table_tes[, xmin := -flank_size]
  peak_table_tes[, xmax := flank_size]
  results_list <- rbind(peak_table_tss, peak_table_tes)
  return(results_list)
}

.score_matrix_RNA <- function(target = NULL, windows = NULL, strand_aware = NULL, weight_col = NULL){
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

.get_translation_boundary_heatmap <- function(index = NULL, query_regions = NULL, tx_length_file = NULL, feature_coords = NULL, flank_size = NULL) {
  sn = names(query_regions[index])
  qregion = query_regions[[index]]
  tx_length_for_translation_site <- data.table::copy(tx_length_file)
  tx_length_for_translation_site <- tx_length_for_translation_site[utr5_len > flank_size & utr3_len > flank_size & cds_len > flank_size, .(tx_name, utr5_len, cds_len, utr3_len)]
  tss_gene_ranges <- tx_length_for_translation_site[, .(seqnames = tx_name, start = utr5_len + 1 - flank_size, end = utr5_len + 1 + flank_size, strand = "*", tx_name = tx_name)]
  tes_gene_ranges <- tx_length_for_translation_site[, .(seqnames = tx_name, start = utr5_len + cds_len - flank_size, end = utr5_len + cds_len + flank_size, strand = "*", tx_name = tx_name)]
  tss_gene_ranges <- GenomicRanges::makeGRangesFromDataFrame(tss_gene_ranges)
  tes_gene_ranges <- GenomicRanges::makeGRangesFromDataFrame(tes_gene_ranges)
  if (length(unique(rtracklayer::width(qregion))) != 1) {
    qregion <- SummarizedExperiment::resize(qregion, 1, fix = "center")
  }
  qregion <- qregion[as.character(GenomicRanges::seqnames(qregion)) %in% as.character(GenomeInfoDb::seqlevels(feature_coords))]
  GenomeInfoDb::seqlevels(qregion) <- GenomicRanges::intersect(GenomeInfoDb::seqlevels(qregion), GenomeInfoDb::seqlevels(feature_coords))
  peak_gene_ranges <- GenomicFeatures::mapToTranscripts(qregion, feature_coords)
  gene_ranges <- list()
  gene_ranges[["tss"]] <- tss_gene_ranges
  gene_ranges[["tes"]] <- tes_gene_ranges
  num_cores <- ifelse(Sys.info()["sysname"] == "Windows", 1, 4)
  sm_all <- suppressWarnings(parallel::mclapply(gene_ranges, .score_matrix_RNA, target = peak_gene_ranges, strand_aware = TRUE, mc.cores = num_cores))
  sm_tss <- sm_all[["tss"]]
  sm_tes <- sm_all[["tes"]]
  return(.peak_heatmap(sm_tss, sm_tes))
}
