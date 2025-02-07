.get_eejunct_new <- function(index = NULL, peak_gr = NULL, tx_by_sec = NULL, flank_size  = NULL) {
  query_regions = peak_gr[[index]]
  testt = BiocGenerics::unlist(tx_by_sec)
  testt = testt[GenomicRanges::width(testt) > flank_size]
  transcript = BiocGenerics::unlist(range(tx_by_sec))

  ss5p <- GenomicRanges::GRanges(
    GenomicRanges::seqnames(testt),
    IRanges::IRanges(ifelse(GenomicRanges::strand(testt) == "+",
                            GenomicRanges::end(testt) + 1,
                            GenomicRanges::start(testt) - 1),
                     ifelse(GenomicRanges::strand(testt) == "+",
                            GenomicRanges::end(testt) + 1,
                            GenomicRanges::start(testt) - 1)),
    GenomicRanges::strand(testt)
  )
  names(ss5p) = names(testt)
  ss3p <- GenomicRanges::GRanges(
    GenomicRanges::seqnames(testt),
    IRanges::IRanges(ifelse(GenomicRanges::strand(testt) == "+",
                            GenomicRanges::start(testt) - 1,
                            GenomicRanges::end(testt) + 1),
                     ifelse(GenomicRanges::strand(testt) == "+",
                            GenomicRanges::start(testt) - 1,
                            GenomicRanges::end(testt) + 1)),
    GenomicRanges::strand(testt)
  )
  names(ss3p) = names(testt)
  transcript_5p <- GenomicRanges::GRanges(
    GenomicRanges::seqnames(transcript),
    IRanges::IRanges(ifelse(GenomicRanges::strand(transcript) == "+",
                            GenomicRanges::end(transcript) + 1,
                            GenomicRanges::start(transcript) -1 ),
                     ifelse(GenomicRanges::strand(transcript) == "+",
                            GenomicRanges::end(transcript) + 1,
                            GenomicRanges::start(transcript) - 1)),
    GenomicRanges::strand(transcript)
  )

  transcript_3p <- GenomicRanges::GRanges(
    GenomicRanges::seqnames(transcript),
    IRanges::IRanges(ifelse(GenomicRanges::strand(transcript) == "+",
                            GenomicRanges::start(transcript) - 1,
                            GenomicRanges::end(transcript) + 1),
                     ifelse(GenomicRanges::strand(transcript) == "+",
                            GenomicRanges::start(transcript) - 1,
                            GenomicRanges::end(transcript) + 1)),
    GenomicRanges::strand(transcript)
  )

  overlaps = GenomicRanges::findOverlaps(ss5p, transcript_5p)
  ss5p = ss5p[-S4Vectors::queryHits(overlaps)]
  overlaps = GenomicRanges::findOverlaps(ss3p, transcript_3p)
  ss3p = ss3p[-S4Vectors::queryHits(overlaps)]


  flanks_ss5p <- GenomicRanges::flank(x = ss5p,
                                      width = flank_size,
                                      start = TRUE,
                                      both = FALSE)
  flanks_ss5p <- GenomicRanges::resize(x = flanks_ss5p,
                                       width = flank_size * 2,
                                       fix = "start")

  flanks_ss3p <- GenomicRanges::flank(x = ss3p,
                                      width = 2 * flank_size,
                                      start = TRUE,
                                      both = FALSE)
  flanks_ss3p <- GenomicRanges::resize(x = flanks_ss3p,
                                       width = flank_size * 3,
                                       fix = "start")

  samplename <- names(peak_gr[index])
  hits_5p <- GenomicFeatures::mapToTranscripts(query_regions, flanks_ss5p, ignore.strand = FALSE)
  hits_5p <- data.table::as.data.table(hits_5p)
  setnames(hits_5p, c('seqnames','xHits'), c('transcriptID','index'))
  peak_table <- data.table::as.data.table(query_regions)
  peak_table[, index:=.I]
  peak_table_5p <- peak_table[hits_5p, on = 'index']
  peak_table_5p[, c("index","transcriptsHits") := NULL]
  peak_table_5p[, Sample := samplename]
  peak_table_5p[, feature := "5PSS"]
  peak_table_5p[, rel_pos := i.start - flank_size]
  peak_table_5p[, flank_size := flank_size]
  hits_3p <- GenomicFeatures::mapToTranscripts(query_regions, flanks_ss3p, ignore.strand = FALSE)
  hits_3p <- data.table::as.data.table(hits_3p)
  setnames(hits_3p, c('seqnames','xHits'), c('transcriptID','index'))
  peak_table_3p <- peak_table[hits_3p, on = 'index']
  peak_table_3p[, c("index", "transcriptsHits") := NULL]
  peak_table_3p[, Sample := samplename]
  peak_table_3p[, feature :="3PSS"]
  peak_table_3p[, rel_pos := i.start - 3 * flank_size]
  peak_table_3p[, flank_size := flank_size]
  return(rbind(peak_table_5p, peak_table_3p))
}
