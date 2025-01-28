.generate_chiped_transcriptome <- function(component)
{
  tx_component_GRange <- GenomicRanges::GRanges()
  for (tx_type in component$tx_types) {
    for (component_type in component[[tx_type]]$component_types) {
      temp_gr <- unlist(component[[tx_type]][[component_type]])
      GenomicRanges::mcols(temp_gr) <- NULL
      GenomicRanges::mcols(temp_gr)$tx_type <- tx_type
      GenomicRanges::mcols(temp_gr)$component_type <- component_type
      tx_component_GRange <- c(tx_component_GRange, temp_gr)
    }
  }
  tx_chiped_GRange <- GenomicRanges::disjoin(tx_component_GRange)

  ret <- list(
    tx_component_GRange = tx_component_GRange,
    tx_chiped_GRange = tx_chiped_GRange
  )

  return(ret)
}

.extract_component <- function(txdb_file = NULL,
                               tx_five_utr_min_length = NULL,
                               tx_cds_min_length = NULL,
                               tx_three_utr_min_length = NULL,
                               tx_long_ncrna_min_length = NULL,
                               tx_promoter_length = NULL,
                               tx_tail_length = NULL,
                               tx_ambiguity = NULL,
                               tx_primary_only = NULL,
                               plot_tx_type =  NULL
)
{
  component <- list()
  component$tx_types <- list()
  tx_lengths <- GenomicFeatures::transcriptLengths(txdb_file)
  print(paste("There are", length(tx_lengths$tx_id), "transcripts of", length(unique(tx_lengths$gene_id)), "genes in the genome."))
  if (tx_primary_only){
    res <- rtracklayer::as.data.frame(tx_lengths %>% group_by(gene_id) %>% filter("tx_len" == max("tx_len")))
    name_filter_tx <- res$tx_name
    print(paste("There are", length(name_filter_tx), "primary transcripts of", length(unique(res$gene_id)), "genes in the genome."))
  }else{
    name_filter_tx <- tx_lengths$tx_name
  }
  tx <- GenomicFeatures::exonsBy(txdb_file, by = "tx", use.names=TRUE)
  print(paste("total", length(tx), "transcripts extracted ..."));

  overlap_count <- GenomicRanges::countOverlaps(tx, tx)
  name_filter_tx <- names(tx[overlap_count < (tx_ambiguity + 2)])
  print(paste("total", length(name_filter_tx), "transcripts left after ambiguity filter ..."))
  tx <- tx[name_filter_tx]
  tx_range <- range(tx)
  name_filter_tx <- names(tx[vapply(tx_range, NROW,numeric(1)) == 1])
  tx <- tx[name_filter_tx]
  print(paste("total", length(name_filter_tx), "transcripts left after check chromosome validity ..."))


  cds <- GenomicFeatures::cdsBy(txdb_file, by = "tx",use.names=TRUE)
  utr5 <- GenomicFeatures::fiveUTRsByTranscript(txdb_file, use.names=TRUE)
  utr3 <- GenomicFeatures::threeUTRsByTranscript(txdb_file, use.names=TRUE)

  utr5_flag <- (sum(width(utr5)) > tx_five_utr_min_length)
  utr5_name <- names(utr5)[utr5_flag]
  utr3_flag <- (sum(width(utr3)) > tx_three_utr_min_length)
  utr3_name <- names(utr3)[utr3_flag]
  cds_flag <- (sum(width(cds)) > tx_cds_min_length)
  cds_name <- names(cds)[cds_flag]
  mRNA_name <- dplyr::intersect(dplyr::intersect(utr5_name,utr3_name),cds_name)
  name_filter_mRNA <- dplyr::intersect(mRNA_name, name_filter_tx)
  print(paste("total",length(name_filter_mRNA),"mRNAs left after component length filter ..."))


  all_mRNA <- IRanges::unique(c(names(utr5),names(utr3),names(cds)))
  ncRNA_name <- dplyr::setdiff(name_filter_tx, all_mRNA)
  ncRNA <- tx[ncRNA_name]
  ncrna_overlap_mrna <- GenomicRanges::countOverlaps(ncRNA, tx[name_filter_mRNA])
  names_overlap_ncrna <- names(ncRNA[ncrna_overlap_mrna < 1])
  ncRNA_flag <-
    (sum(width(ncRNA)) > tx_long_ncrna_min_length)
  names_flag_ncRNA <- names(ncRNA)[ncRNA_flag]
  name_filter_ncRNA <- dplyr::intersect(names_overlap_ncrna, names_flag_ncRNA)
  print(paste("total",length(name_filter_ncRNA),"ncRNAs left after ncRNA length filter ..."))


  tx <- tx[name_filter_tx]
  tx_range <- tx_range[name_filter_tx]


  promoter <- SummarizedExperiment::flank(tx_range, tx_promoter_length, start = TRUE)
  tail <- SummarizedExperiment::flank(tx_range, tx_tail_length, start = FALSE)

  tx_GRange <- BiocGenerics::unlist(tx)
  GenomicRanges::mcols(tx_GRange) <- NULL
  tx_with_flank_gr <- c(BiocGenerics::unlist(promoter), tx_GRange, BiocGenerics::unlist(tail))
  tx_with_flank <- rtracklayer::split(tx_with_flank_gr, names(tx_with_flank_gr))
  tx_with_flank <- GenomicRanges::reduce(tx_with_flank)

  for(tx_types in plot_tx_type )
  {
    if(tx_types == "tx" )
    {
      print("generate components for all tx")
      {
        component$tx_types <- c(component$tx_types, "tx")
        component[[tx_types]]$component_types <- c("promoter", "rna", "tail")
        component[[tx_types]]$names <- name_filter_tx
        component[[tx_types]]$tx_with_flank <- tx_with_flank
        component[[tx_types]]$tx_with_flank_len <- sum(width(component$tx$tx_with_flank))
        component[[tx_types]]$tx <- tx
        component[[tx_types]]$promoter <- promoter
        component[[tx_types]]$rna <- tx
        component[[tx_types]]$tail <- tail
        component[[tx_types]]$tx_range <- tx_range
      }
    }
    if(tx_types == "mRNA" )
    {
      print("generate components for mRNA")
      if (length(name_filter_mRNA) > 0)
      {
        component$tx_types <- c(component$tx_types, "mRNA")
        component[[tx_types]]$component_types <- c("promoter", "utr5", "cds", "utr3", "tail")
        component[[tx_types]]$names <- name_filter_mRNA
        component[[tx_types]]$tx_with_flank <- tx_with_flank[name_filter_mRNA]
        component[[tx_types]]$tx_with_flank_len <- sum(width(component$mRNA$tx_with_flank))
        component[[tx_types]]$tx <- tx[name_filter_mRNA]
        component[[tx_types]]$promoter <- promoter[name_filter_mRNA]
        component[[tx_types]]$utr5 <- utr5[name_filter_mRNA]
        component[[tx_types]]$cds <- cds[name_filter_mRNA]
        component[[tx_types]]$utr3 <- utr3[name_filter_mRNA]
        component[[tx_types]]$tail <- tail[name_filter_mRNA]
      }
    }
    if(tx_types == "ncRNA" )
    {
      print("generate components for lncRNA")
      if (length(name_filter_ncRNA) > 0)
      {
        component$tx_types <- c(component$tx_types, "ncRNA")
        component[[tx_types]]$component_types <- c("promoter", "ncRNA", "tail")
        component[[tx_types]]$names <- name_filter_ncRNA
        component[[tx_types]]$tx_with_flank <- tx_with_flank[name_filter_ncRNA]
        component[[tx_types]]$tx_with_flank_len <- sum(width(component$ncRNA$tx_with_flank))
        component[[tx_types]]$tx <- tx[name_filter_ncRNA]
        component[[tx_types]]$promoter <- promoter[name_filter_ncRNA]
        component[[tx_types]]$ncRNA <- tx[name_filter_ncRNA]
        component[[tx_types]]$tail <- tail[name_filter_ncRNA]
      }
    }
  }
  print("generate chiped transcriptome")
  rslt <- .generate_chiped_transcriptome(component)
  component$tx$tx_component_GRange <- rslt$tx_component_GRange
  component$tx$tx_chiped_GRange <- rslt$tx_chiped_GRange
  return(component)
}

.generate_checking_ranges <- function(tx_information, component_types, checking_ranges_number = 500)
{
  tx_number <- length(tx_information$names)
  component_type_number <- length(component_types)


  component_width <- matrix(0, tx_number, component_type_number)
  rownames(component_width) <- tx_information$names
  colnames(component_width) <- component_types
  end_point <- component_width
  start_point <- component_width

  for (component_type in component_types) {
    component_width[, component_type] <- sum(width(tx_information[[component_type]]))
  }

  component_width_ratio <- component_width / rowSums(component_width)
  component_width_ratio_avg <- colSums(component_width_ratio)
  component_width_average <- floor(component_width_ratio_avg / sum(component_width_ratio_avg) * checking_ranges_number + 0.5)

  end_point <- t(apply(component_width, 1, cumsum))
  start_point[, 1] <- 0
  if (component_type_number > 1) {
    start_point[, seq(2,component_type_number)] <- end_point[, seq_len(component_type_number-1)]
  }
  start_point <- start_point + 1

  component_width_average_mat <- replicate(tx_number, component_width_average)
  if (component_type_number > 1) {
    component_width_average_mat <- t(component_width_average_mat)
  }

  ret <- list(
    component_width = component_width,
    start_point = start_point,
    end_point = end_point,
    component_width_average = component_width_average
  )

  return(ret)
}

.generate_Guitar_coord_txdb <- function(component)
{
  guitar_txdb <- list()
  guitar_txdb$tx_types <- component$tx_types

  for (tx_type in component$tx_types) {
    print(paste("generate coverage checking ranges for", tx_type))


    guitar_txdb[[tx_type]]$tx <- component[[tx_type]]$tx_with_flank
    guitar_txdb[[tx_type]]$tx_length <- sum(width(component[[tx_type]]$tx_with_flank))
    component_types <- component[[tx_type]]$component_types


    component_type_number <- length(component_types)

    rslt <- .generate_checking_ranges(component[[tx_type]], component_types, checking_ranges_number = 500)

    guitar_txdb[[tx_type]]$component_width <- rslt$component_width
    guitar_txdb[[tx_type]]$component_width_ptc <- rslt$component_width/ apply(rslt$component_width, 1, sum)
    guitar_txdb[[tx_type]]$start_point <- rslt$start_point
    guitar_txdb[[tx_type]]$end_point <- rslt$end_point
    guitar_txdb[[tx_type]]$component_width_average <- rslt$component_width_average

    if (tx_type == 'tx') {
      guitar_txdb[[tx_type]]$component_width_average_pct <- guitar_txdb[[tx_type]]$component_width_average / sum(guitar_txdb[[tx_type]]$component_width_average)
    }
    if (tx_type == 'mRNA') {
      guitar_txdb[[tx_type]]$component_width_average_pct <- guitar_txdb[[tx_type]]$component_width_average / sum(guitar_txdb[[tx_type]]$component_width_average)
    }
    if (tx_type == 'ncRNA') {
      guitar_txdb[[tx_type]]$component_width_average_pct <- guitar_txdb[[tx_type]]$component_width_average / sum(guitar_txdb[[tx_type]]$component_width_average)
    }

    cumsum_start_average_pct <- cumsum(guitar_txdb[[tx_type]]$component_width_average_pct)
    guitar_txdb[[tx_type]]$component_start_average_pct[1] <- 0
    if (component_type_number > 1) {
      guitar_txdb[[tx_type]]$component_start_average_pct[seq(2, component_type_number)] <- cumsum_start_average_pct[seq_len(component_type_number - 1)]
      names(guitar_txdb[[tx_type]]$component_start_average_pct) <- component_types
    }
  }

  guitar_txdb$tx$tx_component_GRange <- component$tx$tx_component_GRange
  guitar_txdb$tx$tx_chiped_GRange <- component$tx$tx_chiped_GRange

  return(guitar_txdb)
}

.make_Guitar_txdb <- function(txdb_file = txdb_file,
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
                             plot_tx_type = plot_tx_type
)
{
  tx_component <- .extract_component(txdb_file = txdb_file,
                                     tx_five_utr_min_length = tx_five_utr_min_length,
                                     tx_cds_min_length = tx_cds_min_length,
                                     tx_three_utr_min_length = tx_three_utr_min_length,
                                     tx_long_ncrna_min_length = tx_long_ncrna_min_length,
                                     tx_promoter_length = tx_promoter_length,
                                     tx_tail_length = tx_tail_length,
                                     tx_ambiguity = tx_ambiguity,
                                     tx_primary_only = tx_primary_only,
                                     plot_tx_type = plot_tx_type)

  guitar_txdb <- .generate_Guitar_coord_txdb(tx_component)

  return(guitar_txdb)
}

.get_Guitar_txdb <- function(txdb_file = NULL,
                             tx_five_utr_min_length = NULL,
                             tx_cds_min_length = NULL,
                             tx_three_utr_min_length = NULL,
                             tx_long_ncrna_min_length = NULL,
                             tx_promoter_length = NULL,
                             tx_tail_length = NULL,
                             tx_ambiguity = NULL,
                             txTxComponentProp = NULL,
                             txMrnaComponentProp = NULL,
                             txLncrnaComponentProp = NULL,
                             tx_primary_only = NULL,
                             plot_tx_type = NULL
)
{
  guitar_txdb <- .make_Guitar_txdb(txdb_file = txdb_file,
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
                                  plot_tx_type = plot_tx_type
  )


  return(guitar_txdb)
}

.getStGroup <- function(bed_file = NULL, set_GRange_lists = NULL, set_group_name = NULL
)
{
  if (!(is.null(bed_file))) {
    set_GRange_lists = vector("list", length(bed_file))
    for (i in seq_len(length(bed_file))) {
      print(paste("import BED file", bed_file[[i]], sep = " "))
      set_GRange_lists[[i]] <-  rtracklayer::blocks(rtracklayer::import(bed_file[[i]]))
    }
  }

  set_group_name_length = length(set_GRange_lists)

  if (!(is.null(set_group_name))) {
    set_group_name_length = length(set_group_name)
    names(set_GRange_lists) <- set_group_name
  } else {
    set_group_name <- paste0("Group", seq_len(set_group_name_length))
    set_group_name_length = set_group_name_length
    names(set_GRange_lists) <- set_group_name
  }

  if (set_group_name_length != set_group_name_length) {
    stop(paste0("Site group is  must be assigned by either bed_file or set_GRange_lists"))
  }
  return(set_GRange_lists)
}

.GRanges_list_map_to_transcripts <- function(site = NULL, map_filter_transcript = FALSE, transcripts = NULL)
{
  if (is(site,"CompressedGRangesList"))
  {
    names(site) <- 1:length(site)
    x_widthes <- sum(rtracklayer::width(site))
    names(x_widthes) <- names(site)
    x_unlisted <- unlist(site)
  } else {
    names(site) <- 1:length(site)
    x_widthes <- rtracklayer::width(site)
    names(x_widthes) <- names(site)
    x_unlisted <- site
  }

  tx_coord <- suppressWarnings(GenomicFeatures::mapToTranscripts(x_unlisted, transcripts, ignore.strand = FALSE))

  x_hit_tx_hit_joint <- paste(names(tx_coord), tx_coord$transcriptsHits, sep = '-')
  tx_coord_grouped <- rtracklayer::split(tx_coord, x_hit_tx_hit_joint)
  mapping_reduced <- GenomicRanges::reduce(tx_coord_grouped)

  mapping_reduced_width <- rtracklayer::width(mapping_reduced)
  mapping_region_nums <- lapply(mapping_reduced_width, function(x) length(x))
  index_of_continous <- which(mapping_region_nums == 1)
  mapping_filter <- mapping_reduced[index_of_continous]
  tx_coord_filtered <- unlist(mapping_filter)

  x_hit_tx_hit <- strsplit(names(tx_coord_filtered), '-')
  x_hits <- as.numeric(lapply(x_hit_tx_hit, `[`, 1))
  tx_hits <- as.numeric(lapply(x_hit_tx_hit, `[`, 2))
  GenomicRanges::mcols(tx_coord_filtered) <- data.frame(x_hits = x_hits, tx_hits= tx_hits)

  if(map_filter_transcript) {
    tx_coord_filtered_width <- rtracklayer::width(tx_coord_filtered)
    tx_coord_filtered <- tx_coord_filtered[tx_coord_filtered_width == x_widthes[tx_coord_filtered$x_hits]]
  }
  idx <- GenomicRanges:::get_out_of_bound_index(tx_coord_filtered)
  if (length(idx) != 0) {
    tx_coord_filtered <- tx_coord_filtered[-idx]
  }
  return(tx_coord_filtered)
}

.sample_points <- function(sites_GRange_lists = NULL,
                          set_sample_number = NULL,
                          set_ambiguity = NULL,
                          plot_tx_type = NULL,
                          sample_model = NULL,
                          map_filter_transcript = FALSE,
                          guitar_txdb = NULL)
{
  map_site_GRanges <- list()
  sites_points <- list()
  sites_GRange_data_frame <-list()
  set_sample_number <- 2 * set_sample_number - 1
  for (tx_type in plot_tx_type) {
    map_site_GRanges[[tx_type]] <- .GRanges_list_map_to_transcripts(sites_GRange_lists[[1]], map_filter_transcript, guitar_txdb[[tx_type]]$tx)
    sites_number <- length(map_site_GRanges[[tx_type]])
    sites_width <- width(map_site_GRanges[[tx_type]])

    if(sample_model == "Equidistance")
    {
      my_function <- function(x,i){
        if (i == 1) {
          round(x / 2)
        } else {
          round(seq(1,x - 1,length.out = i))
        }
      }

      sites_points <-t(vapply(sites_width, my_function, i = set_sample_number, numeric(set_sample_number)))
    }
    if(sample_model == "random")
    {
      my_function<-function(x,i){
        if (i == 1) {
          round(x / 2)
        } else {
          a <- AnnotationDbi::sample(x,i,replace = FALSE)
          b <- SummarizedExperiment::sort(a)
        }
      }
      sites_points <- S4Vectors::t(vapply(sites_width, my_function, i = set_sample_number, numeric(set_sample_number)))
    }
    sites_points_vector <- as.numeric(S4Vectors::t(sites_points))
    sites_points_data_frame <- data.frame(chr = rep(SummarizedExperiment::seqnames(map_site_GRanges[[tx_type]]), each = set_sample_number),
                                          start = sites_points_vector + rtracklayer::start(rep(map_site_GRanges[[tx_type]], each = set_sample_number)) - 1, end = sites_points_vector + rtracklayer::start(rep(map_site_GRanges[[tx_type]], each = set_sample_number)))

    sites_GRange_data_frame[[tx_type]] <- GenomicRanges::makeGRangesFromDataFrame(sites_points_data_frame)
    points_overlap_tx <- stats::ave(seq(map_site_GRanges[[tx_type]]), map_site_GRanges[[tx_type]]$x_hits, FUN = length)
    GenomicRanges::mcols(sites_GRange_data_frame[[tx_type]]) <- data.frame(sites_length = c(rep(sites_width,each = set_sample_number)),x_hits =
                                                                             c(rep(map_site_GRanges[[tx_type]]$x_hits,each = set_sample_number)),points_overlap_tx =
                                                                             c(rep(points_overlap_tx,each = set_sample_number)))
  }
  return(sites_GRange_data_frame)
}

.normalize <- function(sites_GRanges = NULL, guitar_txdb = NULL, tx_type = NULL, overlap_index = NULL, site_length_index = NULL)
{
  sites_information <- data.frame()

  sites_points_position_normalize <- list()

  sites_points_position_tx <- list()

  sites_points_position_tx <- rtracklayer::end(sites_GRanges[[tx_type]])

  names(sites_points_position_tx) <- SummarizedExperiment::seqnames(sites_GRanges[[tx_type]])
  #step 1
  start_point_mat <- guitar_txdb[[tx_type]]$start_point[names(sites_points_position_tx), ]
  start_point_differ <- start_point_mat - sites_points_position_tx
  max_which <- function(x)
  {
    max(which(x < 0))
  }
  component_which <- function(x, componet_pct, sites_points_componet)
  {
    sites_points_componet_pct <- componet_pct[x,][sites_points_componet[[x]]]
  }
  sites_points_componet <- apply(start_point_differ, 1, max_which)
  # step 2
  sites_points_position_componet <- start_point_differ[cbind(seq_along(sites_points_componet), sites_points_componet)] * -1
  # step 3
  sites_points_componet_mat <- start_point_mat[cbind(seq_along(sites_points_componet), sites_points_componet)]
  # step 4
  sites_points_componet_width_avg <- guitar_txdb[[tx_type]]$component_width_average_pct[sites_points_componet]
  # step 5
  sites_points_componet_start_pct <- guitar_txdb[[tx_type]]$component_start_average_pct[sites_points_componet]
  # step 6
  component_width_mat <- guitar_txdb[[tx_type]]$component_width[names(sites_points_position_tx), ]
  sites_points_componet_width <- component_width_mat[cbind(seq_along(sites_points_componet), sites_points_componet)]
  # step 7
  sites_points_position_normalize <- sites_points_position_componet / sites_points_componet_width * sites_points_componet_width_avg + sites_points_componet_start_pct
  names(sites_points_position_normalize) <- sites_GRanges[[tx_type]]$x_hits
  #step 8
  sites_componet_pct <- guitar_txdb[[tx_type]]$component_width_ptc[names(sites_points_componet),]
  sites_points_componet_pct <- unlist(lapply( 1:length(sites_points_componet), component_which, sites_componet_pct, sites_points_componet))
  sites_points_weight <-  sites_points_componet_width_avg / (sites_GRanges[[tx_type]]$points_overlap_tx ^ overlap_index)/ sites_points_componet_pct * (sites_GRanges[[tx_type]]$sites_length ^ site_length_index)
  names(sites_points_weight) <- sites_GRanges[[tx_type]]$x_hits
  return(list(sites_points_position_normalize,sites_points_weight))
}

.generate_density_confidence_interval <- function(
    sites_points_relative = NULL,
    sites_point_weight = NULL,
    confidence_interval_resampling_time = NULL,
    adjustment_factor = NULL,
    enable_confidence_interval = enable_confidence_interval,
    confidence_interval_range = c(0.025,0.975))
{
  density_data_frame <- data.frame()
  for (group_name in names(sites_points_relative)) {
    sites_point_weight[[group_name]] <- na.omit(sites_point_weight[[group_name]])
    sites_points_relative[[group_name]] <- na.omit(sites_points_relative[[group_name]])
    sites_weight <- sites_point_weight[[group_name]] / sum(sites_point_weight[[group_name]])
    site_id <- sites_points_relative[[group_name]]
    fit1 <- suppressWarnings(density(site_id, adjustment_factor = adjustment_factor, from = 0, to = 1, n = 256, weight = sites_weight))


    tmp <-  data.frame(
      x = fit1$x,
      density = fit1$y / (sum(diff(fit1$x) * (head(fit1$y, -1) + tail(fit1$y, -1))) / 2),
      group = rep(group_name, times = length(density))
    )

    if (enable_confidence_interval)
    {
      point_ind <- seq(1, length(site_id))
      names(point_ind) <- names(site_id)
      point_ind_grouped <- split(point_ind, names(point_ind))
      fit2 <- replicate(
        confidence_interval_resampling_time,
        {
          point_ind_grouped_resampled <- sample(names(point_ind_grouped), replace = TRUE)
          point_ind_resampled <- unlist(point_ind_grouped[point_ind_grouped_resampled])
          site_sampled_id <- site_id[point_ind_resampled];
          weight_sapmpled <-sites_weight[point_ind_resampled];
          suppressWarnings(density(site_sampled_id, adjustment_factor = adjustment_factor,from = 0, to = 1,n = 256,weight = weight_sapmpled)$y)
        }
      )

      fit3 <- apply(fit2, 1, quantile, confidence_interval_range)

      tmp$confidence_down <- fit3[1, ]
      tmp$confidence_up <- fit3[2, ]
    }

    density_data_frame <- rbind(density_data_frame, tmp)
  }

  return(density_data_frame)
}

.generate_pos_para <- function(peak)
{
  pos <- list()
  pos$fig_top <- 1.05 * peak
  pos$fig_bottom <- -0.1 * peak
  pos$rna_comp_text <- -0.06 * peak
  pos$rna_lgd_bl <- -0.03 * peak
  pos$rna_lgd_h_cds <- 0.01 * peak
  pos$rna_lgd_h_ncrna <- 0.01 * peak
  pos$rna_lgd_h_rna <- 0.01 * peak
  pos$rna_lgd_h_utr <- 0.005 * peak
  pos$rna_lgd_h_flank <- 0.002 * peak

  return(pos)
}

.RNA_plot_structure <- function(p = NULL, component_width = NULL, head_or_tail = NULL, pos = NULL, tx_promoter_length = NULL, tx_tail_length = NULL)
{
  component_structure <- data.frame(width = component_width)
  component_structure$end <- cumsum(component_width)
  component_structure$start <- c(0, component_structure$end[seq_len(length(component_width) - 1)]) + 0.001
  component_structure$mid <- (component_structure$start + component_structure$end) / 2
  component_structure$width <- component_width
  component_structure$comp <- names(component_width)
  component_structure$label <- names(component_width)
  if (head_or_tail)
  {
    component_structure$label[component_structure$comp == "promoter"] <- paste("Promoter(", tx_promoter_length / 1000, "kb)", sep = "")
  } else {
    component_structure$label[component_structure$comp == "promoter"] <- NA
  }
  component_structure$label[component_structure$comp == "utr5"] <- "5'UTR"
  component_structure$label[component_structure$comp == "cds"] <- "CDS"
  component_structure$label[component_structure$comp == "utr5"] <- "5'UTR"
  component_structure$label[component_structure$comp == "ncRNA"] <- "ncRNA"
  component_structure$label[component_structure$comp == "rna"] <- "RNA"
  component_structure$label[component_structure$comp == "utr3"] <- "3'UTR"
  if (head_or_tail)
  {
    component_structure$label[component_structure$comp == "tail"] <-  paste("Tail(", tx_tail_length / 1000, "kb)",sep = "")
  } else {
    component_structure$label[component_structure$comp == "tail"] <-  NA
  }
  component_structure$alpha <- 0.99
  component_structure$alpha[component_structure$comp == "cds"] <- 0.272
  component_structure$alpha[component_structure$comp == "ncRNA"] <- 0.2
  component_structure$alpha[component_structure$comp == "rna"] <- 0.2
  component_structure$lgd_height <- pos$rna_lgd_h_flank
  component_structure$lgd_height[component_structure$comp == "cds"] <- pos$rna_lgd_h_cds
  component_structure$lgd_height[component_structure$comp == "ncRNA"] <- pos$rna_lgd_h_ncrna
  component_structure$lgd_height[component_structure$comp == "rna"] <- pos$rna_lgd_h_rna
  component_structure$lgd_height[component_structure$comp == "utr5"] <- pos$rna_lgd_h_utr
  component_structure$lgd_height[component_structure$comp == "utr3"] <- pos$rna_lgd_h_utr

  for (comp in component_structure$comp) {
    x = component_structure[comp, "mid"]
    y = pos$rna_comp_text
    label = component_structure[comp, "label"]
    p <- p + ggplot2::annotate("text", x = x, y = y, label = label)

    xmin = component_structure[comp, "start"]
    xmax = component_structure[comp, "end"]
    ymin = pos$rna_lgd_bl - component_structure[comp, "lgd_height"]
    ymax = pos$rna_lgd_bl + component_structure[comp, "lgd_height"]
    alpha = component_structure[comp, "alpha"]
    p <- p + ggplot2::annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha, colour = "black")

    xmin = component_structure[comp, "start"]
    xmax = component_structure[comp, "end"]
    ymin = pos$rna_lgd_bl - component_structure[comp, "lgd_height"]
    ymax = pos$rna_lgd_bl + component_structure[comp, "lgd_height"]
    alpha = component_structure[comp, "alpha"]
    p <- p + annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha, colour = "black")
  }

  if (nrow(component_structure) > 1) {
    vline_pos <- data.frame(
      x1 = component_structure[seq_len(nrow(component_structure) -1), "end"],
      x2 = component_structure[seq_len(nrow(component_structure) -1), "end"],
      y1 = rep(pos$fig_top, 4),
      y2 = rep(pos$rna_lgd_bl, 4)
    )
    p <- p + suppressWarnings(
      ggplot2::geom_segment(
        ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
        linetype = "dotted", size = 1, data = vline_pos
      )
    )

  }


  return(p)
}

.plot_density_confidence_interval <- function(density_data_frame_confidence_interval = NULL,
                                              component_width = NULL,
                                              head_or_tail = NULL,
                                              title = NULL,
                                              enable_confidence_interval = NULL,
                                              tx_type = NULL,
                                              tx_promoter_length = NULL,
                                              tx_tail_length = NULL)
{
  if (enable_confidence_interval) {
    peak <- max(density_data_frame_confidence_interval$confidence_up)
  } else {
    peak <- max(density_data_frame_confidence_interval$density)
  }
  pos <- .generate_pos_para(peak)

  samples <- factor(density_data_frame_confidence_interval$group)
  p <- suppressWarnings(ggplot2::ggplot(density_data_frame_confidence_interval, ggplot2::aes(x = x)))
  p <- p + ggplot2::geom_line(ggplot2::aes(y = density, colour = samples), alpha = 1, size = 1)
  p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = rep(0, length(density)), ymax = density, colour = samples, fill = samples), alpha = 0.2 )
  if (enable_confidence_interval) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = density_data_frame_confidence_interval$confidence_down, colour = samples), colour="blue", alpha = 0.4, size = 0.5)
    p <- p + ggplot2::geom_line(ggplot2::aes(y = density_data_frame_confidence_interval$confidence_up, colour = samples), colour="black", alpha = 0.4, size = 0.5)
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = density_data_frame_confidence_interval$confidence_down, ymax = density_data_frame_confidence_interval$confidence_up, colour = samples, fill = samples), alpha = 0.2, colour = NA)
  }

  if (tx_type == "tx"){
    p <- p + ggplot2::labs(title = paste("Distribution on exon")) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                        margin = margin(10, 0, 10, 0),
                                                        size = 14))
  }
  if (tx_type == "mRNA"){
    p <- p + ggplot2::labs(title = paste("Distribution on mRNA")) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                        margin = margin(10, 0, 10, 0),
                                                        size = 14))
  }
  if (tx_type == "ncRNA"){
    p <- p + ggplot2::labs(title = paste("Distribution on ncRNA")) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                        margin = ggplot2::margin(10, 0, 10, 0),
                                                        size = 14))
  }
  p <- p + ggplot2::theme(axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                          axis.title.x = ggplot2::element_blank())
  p <- p + ggplot2::ylab("Frequency")
  p <- p + ggplot2::theme(axis.title.y = ggplot2::element_text(vjust = 2, size = 14,face = "bold"))
  p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank(),axis.text.y = ggplot2::element_text(size = 13))
  p <- p + ggplot2::theme(panel.grid.minor = element_blank(),panel.grid.major = ggplot2::element_line(color = "white", size = .7))
  p <- p + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "#4c4c4c", fill = NA, size = 0.8))


  p <- .RNA_plot_structure(p = p, component_width = component_width, head_or_tail = head_or_tail, pos = pos, tx_promoter_length = tx_promoter_length, tx_tail_length = tx_tail_length)

  p <- p + ggplot2::theme(legend.position = "bottom") + ggplot2::theme(legend.text = ggplot2::element_text(size = 13))
  p <- p + ggplot2::theme(legend.title = ggplot2::element_blank())
  p <- p + ggplot2::scale_y_continuous(limits = c(pos$fig_bottom, pos$fig_top),
                                       expand = c(0, 0))

  return(p)
}

.bed_file_test <- function(bed_file) {
  if (tools::file_ext(bed_file) != "bed") {
    stop("The file is not in .bed format")
  }
  first_line <- readLines(bed_file, n = 1)
  if (any(!grepl("^\\d+$", unlist(strsplit(first_line, "\t"))))) {
    warning("This file contains column names. We have removed the column names. Please use a BED file without column names.")
    bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE, skip = 1)
  } else {
    bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
  }
  if (any(!grepl("^\\d+$", bed_data$V2))) {
    non_numeric_rows <- which(!grepl("^\\d+$", bed_data$V2))
    stop(paste("The second column has non-numeric values, error at row(s): ", paste(non_numeric_rows, collapse = ", ")))
  }
  if (any(!grepl("^\\d+$", bed_data$V3))) {
    non_numeric_rows <- which(!grepl("^\\d+$", bed_data$V3))
    stop(paste("The third column has non-numeric values, error at row(s): ", paste(non_numeric_rows, collapse = ", ")))
  }
  if (any(bed_data$V3 <= bed_data$V2)) {
    invalid_rows <- which(bed_data$V3 <= bed_data$V2)
    stop(paste("The third column is not greater than the second column, error at row(s): ", paste(invalid_rows, collapse = ", ")))
  }

  return("BED file test passed")
}




