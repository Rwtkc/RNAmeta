#' @title Get Preliminary Analysis Data
#' @param bed_file A list of BED file paths or GRanges objects. Each entry in the list
#' represents a set of genomic regions for analysis.
#' @param txdb_file A file path to a TxDb (transcript database) file. The TxDb file
#' contains genomic annotations related to transcripts, exons, and other genomic features.
#' @param tx_length_file A file path to a data file containing transcript length data.
#' This file should include information about the lengths of different transcripts.
#' @param gff_file A file path to a GFF (General Feature Format) file that provides
#' genomic annotations, including gene and transcript information. The file should
#' contain columns such as gene ID and feature types (e.g., "gene", "CDS", "exon").
#'
#' @returns A list containing the following elements:
#' - `cct_list`: A list of results from batch genomic annotation processes for each
#' BED file.
#' - `peak_gr_list`: A list of GenomicRanges objects representing peaks for each group.
#' - `gff_file`: The loaded GFF file, processed and adjusted for further analysis.
#' - `peak_gr_for_motifscan`: A list of GenomicRanges objects resized for motif scanning.
#' - `motif_parameter`: A vector containing the motif parameter extracted from the input.
#' - `flanking_size`: The flanking size (set to 1000 by default) used for genomic analyses.
#' - `txdb_features`: A list of features extracted from the TxDb file, such as exons,
#' promoters, and coding sequences.
#' - `gene_max_tx_lengths`: A data table containing the maximum transcript lengths
#' for each gene.
#' - `short_flanking_size`: The short flanking size (set to 100 by default) used in
#' certain genomic analysis steps.
#'
#' @export
#'
get_preliminary_analysis_data <- function(bed_file = NULL, txdb_file = NULL, tx_length_file = NULL, gff_file = NULL){
  number_of_bins <- 100
  genomic_annotation_priority <- c("stop_codon", "start_codon", "CDS", "5UTR", "3UTR", "Intron", "Transcript", "Promoter", "Intergenic")
  motif_parameter <- "6,7,8_1000"
  motif_parameter <- BiocGenerics::unlist(strsplit(motif_parameter, "_"))
  subset_for_motif <- motif_parameter[2]
  cct_list <- list()
  peak_gr_list <- list()
  flanking_size <- 1000
  short_flanking_size <- 100
  union_intersection <- "other"
  bed_file_list <- list()
  bed_file_list <- bed_file
  txdb_file_path = txdb_file
  txdb_file <- AnnotationDbi::loadDb(txdb_file)
  tx_length_file <- load(tx_length_file)
  assign("tx_length_file", get(tx_length_file))
  gene_max_tx_lengths <- tx_length_file[, IRanges::lapply(.SD, max), by = gene_id]
  txdb_features <- .get_txdb_features_from_GRanges(txdb_file = txdb_file, tx_length_file = tx_length_file, number_of_bins = number_of_bins)
  gff_file <- load(gff_file)
  assign("gff_file", get(gff_file))
  if ("gene_biotype" %in% names(GenomicRanges::mcols(gff_file))) {
    names(GenomicRanges::mcols(gff_file))[which(names(GenomicRanges::mcols(gff_file)) == "gene_biotype")] <- "gene_type"
  }
  if ("transcript_biotype" %in% names(GenomicRanges::mcols(gff_file))) {
    names(GenomicRanges::mcols(gff_file))[which(names(GenomicRanges::mcols(gff_file)) == "transcript_biotype")] <- "transcript_type"
  }
  bed_file_lenght <- length(bed_file)
  for(i in 1:bed_file_lenght){
    cct <- NULL
    peak_dt = list()
    peak_gr = list()
    peak_gr_for_motifscan = list()
    group_name <- paste0("Group",i)
    ctt <- .batch_peak_reader(group_name = group_name, subset_for_motif = subset_for_motif, txdb_file = txdb_file, bed_file = bed_file_list[[i]], txdb_file_path = txdb_file_path)
    if(max(GenomicRanges::width(ctt)) < 15) {
      peak_gr_for_motifscan[[group_name]] <- GenomicRanges::resize(ctt, width = 15, fix="center")
    } else {
      peak_gr_for_motifscan[[group_name]] <- ctt
    }
    ctt <- GenomicRanges::resize(ctt, width = 1, fix = "center")
    peak_gr[[group_name]] = ctt
    peak_dt[[group_name]] = data.table::as.data.table(ctt)
    cct = lapply(1:length(peak_gr),
                 .batch_genomic_annotation,
                 peak_gr = peak_gr,
                 txdb_features = txdb_features,
                 genomic_annotation_priority = genomic_annotation_priority,
                 tx_length_file = tx_length_file,
                 gff_file = gff_file,
                 peak_dt = peak_dt,
                 group_name = group_name)
    cct_list <- append(cct_list, cct)
    peak_gr_list <- append(peak_gr_list, peak_gr)
  }
  preliminary_analysis_data <- list(cct_list,
                                    peak_gr_list,
                                    gff_file,
                                    peak_gr_for_motifscan,
                                    motif_parameter,
                                    flanking_size,
                                    txdb_features,
                                    gene_max_tx_lengths,
                                    short_flanking_size)
  return(preliminary_analysis_data)
}










