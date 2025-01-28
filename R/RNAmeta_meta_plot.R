#' @title Plot genomic feature distributions based on sampled points.
#'
#' @param bed_file Path to the BED file defining genomic regions (sites).
#'   Requires standard BED format columns: chrom, start, end, and optional metadata fields.
#' @param txdb_file TxDb object or file path containing transcript annotations.
#'   Supports Bioconductor TxDb packages or SQLite database files with genomic
#'   features including gene IDs, transcript coordinates, and exon structures.
#' @param set_group_name Group names for the genomic sites.
#'   Specifies the names used to categorize or group different genomic regions
#'   defined in the BED file. This parameter allows users to assign meaningful
#'   labels to different sets of sites based on their biological or experimental
#'   significance.
#' @param set_ambiguity Maximum allowed overlap between genomic sites.
#'   Specifies the maximum number of base pairs that can overlap between
#'   different genomic regions defined in the BED file. This parameter helps
#'   control the level of ambiguity when determining overlap between sites.
#'   Default: 5
#' @param set_sample_number The number of bases sampled at each genomic site.
#'   Specifies how many bases are selected for analysis or sampling at each
#'   site defined in the BED file. This parameter controls the granularity
#'   of sampling across the genomic regions.
#'   Default: 10
#' @param sample_model The sampling model to be used.
#'   Specifies the method for sampling bases from each genomic site.
#'   Options include "Equidistance" for evenly spaced sampling across the
#'   genomic region, or "random" for randomly selecting bases within the region.
#'   Default: "Equidistance"
#' @param tx_five_utr_min_length Minimum length of the 5' UTR.
#'   Specifies the minimum length for the 5' untranslated region (UTR) of
#'   transcripts. This parameter helps filter out transcripts with shorter
#'   5' UTRs, allowing for a more focused analysis of genes with longer
#'   regulatory regions.
#'   Default: 100
#' @param tx_cds_min_length Minimum length of the CDS (coding sequence).
#'   Specifies the minimum length for the coding sequence (CDS) of transcripts.
#'   This parameter allows filtering of transcripts with shorter CDS, ensuring
#'   the analysis focuses on genes with more substantial coding regions.
#'   Default: 100
#' @param tx_three_utr_min_length Minimum length of the 3' UTR.
#'   Specifies the minimum length for the 3' untranslated region (UTR) of
#'   transcripts. This parameter helps filter out transcripts with shorter
#'   3' UTRs, allowing for a more focused analysis of genes with longer
#'   regulatory regions.
#'   Default: 100
#' @param tx_long_ncrna_min_length Minimum length of long non-coding RNA (lncRNA).
#'   Specifies the minimum length for long non-coding RNAs (lncRNAs). This
#'   parameter allows filtering of shorter lncRNAs, focusing the analysis
#'   on longer non-coding transcripts that may have more significant biological roles.
#'   Default: 100
#' @param tx_lncrna_overlap_mrna Whether to allow overlap between lncRNA and mRNA.
#'   Specifies whether long non-coding RNAs (lncRNAs) are allowed to overlap
#'   with messenger RNAs (mRNAs) in the analysis. If set to TRUE, lncRNAs and
#'   mRNAs can overlap in their genomic regions. If set to FALSE, any overlap
#'   between lncRNA and mRNA regions will be excluded from the analysis.
#'   Default: FALSE
#' @param tx_promoter_length Length of the promoter region.
#'   Specifies the length of the promoter region to be considered in the
#'   analysis. This parameter defines the region upstream of the transcription
#'   start site (TSS) that is included for promoter analysis.
#'   Default: 1000
#' @param tx_tail_length Length of the tail region.
#'   Specifies the length of the tail region to be considered in the
#'   analysis. This parameter defines the region downstream of the transcription
#'   end site (TES) that is included for tail analysis.
#'   Default: 1000
#' @param tx_ambiguity Maximum overlap between transcripts.
#'   Specifies the maximum number of base pairs that can overlap between
#'   different transcripts. This parameter helps control the level of ambiguity
#'   when determining overlap between transcript regions.
#'   Default: 5
#' @param tx_primary_only Whether to use only the primary transcript.
#'   If set to TRUE, only the primary transcript will be considered for analysis.
#'   If set to FALSE, all transcripts associated with a gene will be included.
#'   Default: FALSE
#' @param map_filter_transcript Whether to filter the length of transcripts to match the original site.
#'   If set to TRUE, only transcripts whose lengths match the length of the original site will be kept.
#'   Default: TRUE
#' @param head_or_tail Whether to retain the promoter (head) and tail regions.
#'   If set to TRUE, both the promoter and tail regions will be considered in the analysis.
#'   If set to FALSE, these regions will be excluded.
#'   Default: TRUE
#' @param enable_confidence_interval Whether to add a confidence interval (CI) curve.
#'   If set to TRUE, a CI curve will be included in the plot.
#'   If set to FALSE, no CI curve will be drawn.
#'   Default: TRUE
#' @param plot_tx_type Which transcript type to plot.
#'   Specifies which transcript type (e.g., "tx", "mRNA", "ncRNA") should be drawn on the plot.
#'   If the specified transcript type is not found in the genome, it will not be plotted.
#'   Default: c("tx", "mRNA", "ncRNA")
#' @param overlap_index Index of site overlap count.
#'   Specifies the number of overlaps between genomic sites to be considered.
#'   Default: 1
#' @param site_length_index Index of site length.
#'   Specifies the index used to filter or categorize genomic sites based on their lengths.
#'   Default: 1
#' @param adjustment_factor Smoothing level of the curve.
#'   Controls the smoothness of the curve drawn in the analysis. Higher values result in a smoother curve.
#'   Default: 1
#' @param confidence_interval_resampling_time Resampling times in density drawing mode.
#'   Specifies the number of resampling iterations to perform when drawing the confidence interval curve.
#'   Default: 1000
#' @param confidence_interval_range Range for the confidence interval.
#'   Defines the upper and lower bounds of the confidence interval for the plot.
#'   Default: c(0.025, 0.975)
#'
#' @returns
#'   Returns a plot object showing the distribution of genomic features (such as mRNA, ncRNA, or general transcripts)
#'   based on the sampled points from the provided BED file and transcript annotations.
#'   The plot includes an optional confidence interval curve if specified.
#'   The plot type and content vary depending on the 'plot_tx_type' parameter and the calculated values
#'   for each genomic region sampled from the BED file.
#' @export
#'
RNAmeta_meta_plot <- function(bed_file = NULL,
                              txdb_file = NULL,
                              set_group_name = NULL,
                              set_ambiguity = 5,
                              set_sample_number = 10,
                              sample_model = c("Equidistance", "random"),
                              tx_five_utr_min_length = 100,
                              tx_cds_min_length = 100,
                              tx_three_utr_min_length = 100,
                              tx_long_ncrna_min_length = 100,
                              tx_lncrna_overlap_mrna = TRUE,
                              tx_promoter_length = 1000,
                              tx_tail_length = 1000,
                              tx_ambiguity = 5,
                              tx_primary_only = TRUE,
                              map_filter_transcript = TRUE,
                              head_or_tail = TRUE,
                              enable_confidence_interval = FALSE,
                              plot_tx_type = c("mRNA", "ncRNA", "tx"),
                              overlap_index = 1,
                              site_length_index = 1,
                              adjustment_factor = 1,
                              confidence_interval_resampling_time = 1000,
                              confidence_interval_range = c(0.025, 0.975)
                              ){
  txdb_file <- AnnotationDbi::loadDb(txdb_file)
  guitar_txdb <- .get_Guitar_txdb(txdb_file = txdb_file,
                                 tx_five_utr_min_length = tx_five_utr_min_length,
                                 tx_cds_min_length = tx_cds_min_length,
                                 tx_three_utr_min_length = tx_three_utr_min_length,
                                 tx_long_ncrna_min_length = tx_long_ncrna_min_length,
                                 tx_promoter_length = tx_promoter_length,
                                 tx_tail_length = tx_tail_length,
                                 tx_ambiguity = tx_ambiguity,
                                 txTxComponentProp = NULL,
                                 txMrnaComponentProp = NULL,
                                 txLncrnaComponentProp = NULL,
                                 tx_primary_only = tx_primary_only,
                                 plot_tx_type = plot_tx_type)
  sites_group <- .getStGroup(bed_file = bed_file, set_GRange_lists = NULL, set_group_name = set_group_name)
  group_names <- names(sites_group)
  sites_group_number <- length(sites_group)
  sites_points_normlize <- list()
  sites_points_relative <- list()
  point_weight <- list()
  for (i in seq_len(sites_group_number)) {
    group_name = group_names[[i]]
    print(paste("sample", set_sample_number, "points for" , group_name, sep = " "))
    sites_points <- sample_points(sites_group[i],
                                  set_sample_number = set_sample_number,
                                  set_ambiguity = set_ambiguity,
                                  plot_tx_type =  plot_tx_type,
                                  sample_model = sample_model,
                                  map_filter_transcript = map_filter_transcript,
                                  guitar_txdb)

    for (tx_type in plot_tx_type) {
      sites_points_normlize[[tx_type]][[group_name]] <- normalize(sites_points, guitar_txdb, tx_type, overlap_index, site_length_index)
      sites_points_relative[[tx_type]][[group_name]] <- sites_points_normlize[[tx_type]][[group_name]][[1]]
      point_weight[[tx_type]][[group_name]] <- sites_points_normlize[[tx_type]][[group_name]][[2]]
    }
  }
  for (tx_type in plot_tx_type)
  {
    print(paste("start figure plotting for", tx_type, "..."))
    if (tx_type == "mRNA") {
      tx_type_name = "mRNA"
    } else if (tx_type == "ncRNA") {
      tx_type_name = "ncRNA"
    } else {
      tx_type_name = "Transcript"
    }
    title <- paste("Distribution on", tx_type_name)
    density_data_frame_confidence_interval <- .generate_density_confidence_interval(
                                                                  sites_points_relative = sites_points_relative[[tx_type]],
                                                                  sites_point_weight = point_weight[[tx_type]],
                                                                  confidence_interval_resampling_time = confidence_interval_resampling_time,
                                                                  adjustment_factor = adjustment_factor,
                                                                  enable_confidence_interval = enable_confidence_interval,
                                                                  confidence_interval_range = confidence_interval_range)
    p <- .plot_density_confidence_interval(
                         density_data_frame_confidence_interval = density_data_frame_confidence_interval,
                         component_width = guitar_txdb[[tx_type]]$component_width_average_pct,
                         head_or_tail = head_or_tail,
                         title = title,
                         enable_confidence_interval = enable_confidence_interval,
                         tx_type = tx_type,
                         tx_promoter_length = tx_promoter_length,
                         tx_tail_length = tx_tail_length)
    return(p)
}}
