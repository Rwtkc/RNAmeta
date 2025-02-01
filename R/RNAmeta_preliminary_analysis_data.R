preliminary_analysis_data <- function(){
  number_of_bins <- 100
  genomic_annotation_priority <- c("Stopcodon", "Startcodon", "CDS", "5UTR", "3UTR", "Intron", "Transcript", "Promoter", "Intergenic")
  motif_parameter <- "6,7,8_1000"
  motif_parameter <- unlist(strsplit(motif_parameter, "_"))
  subset_for_motif_value <- motif_parameter[2]
  cct_list <- list()
  peak_gr_list <- list()
  flanking_size <- 1000
  short_flanking_size <- 100
  union_intersection <- "other"

  bed_file <- list()
  bedfile2 <- "/public/shiny/RNAmeta/data/exampleDataFile/ac_diff_down.bed6.bed"
  txdb <- loadDb("/public/shiny/RNAmeta/data/txdb/ara_TAIR10.txdb.sqlite")
  load("/public/shiny/RNAmeta/data/txlens/ara_TAIR10.txlens.rda")
  txlens2 <- txlens[, lapply(.SD, max), by = gene_id]
  txdbFeatures <- getTxdbFeaturesFromGRanges(txdb, txlens2, numOfBins)
  load("/public/shiny/RNAmeta/data/gff/ara_TAIR10.gff.rda")
  if ("gene_biotype" %in% names(mcols(gff))) {
    names(mcols(gff))[which(names(mcols(gff)) == "gene_biotype")] <- "gene_type"
  }
  if ("transcript_biotype" %in% names(mcols(gff))) {
    names(mcols(gff))[which(names(mcols(gff)) == "transcript_biotype")] <- "transcript_type"
  }
  species <- "ara_TAIR10"


}
