.peak_type <- function(query_regions = NULL, gff_file = NULL) {
  overlaps = GenomicRanges::findOverlaps(query_regions, gff_file)
  overlaps_query = query_regions[S4Vectors::queryHits(overlaps)]
  overlaps_gff = gff_file[S4Vectors::subjectHits(overlaps)]
  overlaps_gff$overlapping_query = paste(GenomicRanges::seqnames(overlaps_query),
                                         GenomicRanges::start(overlaps_query),
                                         GenomicRanges::end(overlaps_query),
                                         GenomicRanges::strand(overlaps_query),
                                         sep=':')
  return(overlaps_gff)
}

.batch_peak_type <- function(index  = NULL, peak_gr  = NULL, gff_file = NULL) {
  peak_gene_types = .peak_type(peak_gr[[index]], gff_file)
  gene_type_stat = data.table::as.data.table(IRanges::table(peak_gene_types$transcript_type))
  gene_type_stat[, Sample := names(peak_gr[index])]
  return(gene_type_stat)
}
