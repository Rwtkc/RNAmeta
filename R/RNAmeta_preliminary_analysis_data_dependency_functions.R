.get_txdb_features_from_GRanges <- function (txdb_file = NULL, tx_length_file = NULL, number_of_bins = NULL) {
  transcripts <- GenomicFeatures::transcripts(txdb_file)
  exons <- GenomicFeatures::exonsBy(x = txdb_file, by = "tx", use.names = TRUE)
  introns <- GenomicFeatures::intronsByTranscript(x = txdb_file, use.names = TRUE)
  promoters <- GenomicFeatures::promoters(txdb_file, upstream = 1000, downstream = 0)
  five_utrs <- GenomicFeatures::fiveUTRsByTranscript(x = txdb_file, use.names = TRUE)
  up_stream = GenomicRanges::flank(five_utrs, 500, start = TRUE, both = FALSE)
  three_utrs <- GenomicFeatures::threeUTRsByTranscript(x = txdb_file, use.names = TRUE)
  down_stream = GenomicRanges::flank(three_utrs, 500, start = FALSE, both = FALSE)
  cds <- GenomicFeatures::cdsBy(x = txdb_file, by = "tx", use.names = TRUE)
  cds_list = range(cds)
  stop_codon = GenomicRanges::flank(cds_list, 200, start = FALSE, both = TRUE)
  start_codon = GenomicRanges::flank(cds_list, 200, start = TRUE, both = TRUE)
  flag_utr5 <- (sum(GenomicRanges::width(five_utrs)) > number_of_bins)
  name_utr5 <- names(five_utrs)[flag_utr5]
  flag_utr3 <- (sum(GenomicRanges::width(three_utrs)) > number_of_bins)
  name_utr3 <- names(three_utrs)[flag_utr3]
  flag_cds <- (sum(GenomicRanges::width(cds)) > number_of_bins)
  name_cds <- names(cds)[flag_cds]
  name_mRNA <- intersect(intersect(name_utr5, name_utr3), name_cds)
  name_filtered_mRNA <- intersect(name_mRNA, names(exons))
  length_filtered_mRNA <- tx_length_file$tx_name
  txdb_features = list(
    'transcripts' = transcripts,
    'exons'       = exons,
    'promoters'   = promoters,
    'five_utrs'    = five_utrs,
    'up_stream'    = up_stream,
    'introns'     = introns,
    'cds'         = cds,
    'three_utrs'   = three_utrs,
    'stop_codon'   = stop_codon,
    'start_codon'  = start_codon,
    'down_stream'  = down_stream
  )
  return(txdb_features)
}

.bed_file_test_for_base <- function(bed_file = NULL, chrom_info = NULL) {
  for (i in 1:length(bed_file)) {
    current_bed_file <- bed_file[[i]]
    if (tools::file_ext(current_bed_file) != "bed") {
      stop(paste("The file", current_bed_file, "is not in .bed format"))
    }
    first_line <- readLines(current_bed_file, n = 1)
    first_five_lines <- readLines(current_bed_file, n = 5)
    possible_delimiters <- c("\t", " ", ",")
    detected_delimiter <- NULL
    for (delimiter in possible_delimiters) {
      split_counts <- sapply(first_five_lines, function(line) length(unlist(strsplit(line, delimiter))))
      if (length(unique(split_counts)) == 1) {
        detected_delimiter <- delimiter
        break
      }
    }
    if (is.null(detected_delimiter)) {
      stop(paste("Unable to detect the delimiter in the BED file", current_bed_file, ". Please check the file format."))
    }
    first_line <- first_five_lines[1]
    first_line <- unlist(strsplit(first_line, detected_delimiter))
    is_header <- suppressWarnings(all(is.na(as.numeric(first_line))))
    if (is_header) {
      warning(paste("This file", current_bed_file, "contains column names. We have removed the column names. Please use a BED file without column names."))
      bed_data <- read.table(current_bed_file, header = FALSE, stringsAsFactors = FALSE, skip = 1)
    } else {
      bed_data <- read.table(current_bed_file, header = FALSE, stringsAsFactors = FALSE)
    }
    if (any(!grepl("^\\d+$", bed_data$V2))) {
      non_numeric_rows <- which(!grepl("^\\d+$", bed_data$V2))
      stop(paste("The second column has non-numeric values in", current_bed_file, "error at row(s):", paste(non_numeric_rows, collapse = ", ")))
    }
    if (any(!grepl("^\\d+$", bed_data$V3))) {
      non_numeric_rows <- which(!grepl("^\\d+$", bed_data$V3))
      stop(paste("The third column has non-numeric values in", current_bed_file, "error at row(s):", paste(non_numeric_rows, collapse = ", ")))
    }
    if (any(bed_data$V3 <= bed_data$V2)) {
      invalid_rows <- which(bed_data$V3 <= bed_data$V2)
      stop(paste("The third column is not greater than the second column in", current_bed_file, "error at row(s):", paste(invalid_rows, collapse = ", ")))
    }
    bed_chroms <- unique(bed_data$V1)
    chrom_info_chroms <- unique(chrom_info$chrom)
    missing_chroms <- setdiff(bed_chroms, chrom_info_chroms)
    if (length(missing_chroms) > 0) {
      stop(paste("The following chromosomes in", current_bed_file, "are not present in chrom_info:", paste(missing_chroms, collapse = ", ")))
    }
    bed_data <- data.table::as.data.table(bed_data)
    if (ncol(bed_data) == 3) {
      bed_data[, V4 := paste0("site", .I)]
      bed_data[, V5 := "."]
      bed_data[, V6 := "*"]
    } else if (ncol(bed_data) == 4) {
      bed_data[, V5 := "."]
      bed_data[, V6 := "*"]
    } else if (ncol(bed_data) == 5) {
      bed_data[, V6 := "*"]
    }
    suppressWarnings(if (sum(!is.na(as.numeric(bed_data[[5]]))) > sum(sapply(bed_data[, 5, with = FALSE], function(col) sum(is.na(as.numeric(col)))))) {
      bed_data <- bed_data[!is.na(as.numeric(bed_data[[5]])), ]
    })
    bed_data <- bed_data[, 1:6, with = FALSE]
  }

  return(bed_data)
}

.txdb_file_test <- function(txdb_file = NULL){
  print("Checking txdb_file begins.")
  if (!grepl("\\.sqlite$", txdb_file)) {
    stop("The file extension of txdb_file must be '.sqlite'")
  }
  txdb_sql <- DBI::dbConnect(RSQLite::SQLite(), dbname = txdb_file)
  if (!"chrominfo" %in% DBI::dbListTables(txdb_sql)) {
    DBI::dbDisconnect(txdb_sql)
    stop("This database file does not meet the query specifications: 'chrominfo' table is missing.")
  }
  chrom_info <- data.table::as.data.table(DBI::dbGetQuery(txdb_sql, "SELECT * FROM chrominfo"))
  if(is.null(chrom_info)){
    DBI::dbDisconnect(txdb_sql)
    stop("This database file does not meet the query specifications: 'chrominfo' table is null.")
  }
  DBI::dbDisconnect(txdb_sql)
  print("txdb_file check passed.")
  return(chrom_info)
}

.batch_peak_reader <- function(group_name = NULL, subset_for_motif = NULL, txdb_file = NULL, bed_file = NULL, txdb_file_path = NULL) {
  group_name = names(group_name)
  bed_headers = c("chr","start","end","name","score","strand")
  peak_file = bed_file
  chrom_info <-.txdb_file_test(txdb_file = txdb_file_path)
  peak_df <- .bed_file_test_for_base(bed_file = peak_file, chrom_info = chrom_info)
  data.table::setnames(peak_df, bed_headers[1:ncol(peak_df)])
  peak_df[!strand %in% c("+", "-", "*"), strand := "*"]
  peak_df[,start := start + 1]
  peak_df <- peak_df[grep("chr", chr, ignore.case = TRUE)]
  peak_df <- peak_df[,chr := gsub("chr", 'chr', chr, ignore.case = TRUE)]
  peak_df <- peak_df[chr %in% GenomeInfoDb::seqlevels(txdb_file)]
  peak_rg = GenomicRanges::GRanges(seqnames = peak_df[,chr], ranges = IRanges::IRanges(peak_df[,start], peak_df[,end]))
  if(ncol(peak_df) >= 6){
    BiocGenerics::strand(peak_rg) <- peak_df$strand
  }
  cn <- colnames(peak_df)
  if (length(cn) > 3) {
    for (j in 4:length(cn)) {
      if(!cn[j] %in% c("strand")) {
        GenomicRanges::mcols(peak_rg)[[cn[j]]] <- peak_df[, S4Vectors::eval(as.name(cn[j]))]
      }
    }
  }
  peak_rg$group_name <- group_name
  peak_rg$peak_name <- paste0("Locus_",1:nrow(peak_df))
  return(IRanges::unique(peak_rg))
}

.get_first_hit_index <- function(x = NULL) {
  which(!duplicated(x))
}

.get_genomic_annotation_internal <- function(peaks = NULL, genomic_region = NULL, type = NULL, same_strand = FALSE){

  if(class(genomic_region) == "GRangesList" | class(genomic_region) == "CompressedGRangesList"){
    GRegion <- BiocGenerics::unlist(genomic_region)
    GRegion$tx_name <- names(GRegion)
    GRegion$length = rtracklayer::width(GRegion)
  } else {
    GRegion <- genomic_region
    GRegion$length = rtracklayer::width(GRegion)
  }

  if (type == "Intron" || type =="CDS") {
    gr2 <- GRegion[!SummarizedExperiment::duplicated(GRegion$tx_name)]
    temp = data.table::data.table(name = GRegion$tx_name)
    temp = temp[,.N,by = name]

    strd <- as.character(rtracklayer::strand(gr2))
    GRegion_length <- temp[,N]
    names(GRegion_length) <- temp[,name]

    GRegion$rank <- lapply(seq_along(strd), function(i) {
      rank <- seq(1, GRegion_length[i])
      if (strd[i] == '-')
        rank <- rev(rank)
      return(rank)
    }) %>% unlist
  }

  if (same_strand) {
    GRegion_hit <- SummarizedExperiment::findOverlaps(peaks, GRegion)
  } else {
    GRegion_hit <- SummarizedExperiment::findOverlaps(peaks, BiocGenerics::unstrand(GRegion))
  }

  if (length(GRegion_hit) == 0) {
    return(NA)
  }

  qh <- S4Vectors::queryHits(GRegion_hit)
  hit.idx <- .get_first_hit_index(x = qh)
  GRegion_hit <- GRegion_hit[hit.idx]
  query_index <- S4Vectors::queryHits(GRegion_hit)
  subject_index <- S4Vectors::subjectHits(GRegion_hit)

  hits <- GRegion[subject_index]
  gene_id <- hits$tx_name

  if (type == "Intron") {
    anno <- paste(type, " (", gene_id, ", intron ", hits$rank,  " of ", GRegion_length[gene_id], ")", sep="")
  } else if (type == "CDS") {
    anno <- paste(type, " (", gene_id, ", CDS ", hits$rank, " of ", GRegion_length[gene_id],", ", hits$length, "bp)", sep="")
  } else if (type == "five_utr") {
    anno <- paste0("UTR5", " (", gene_id, ")")
  } else if (type == "three_utr") {
    anno <- paste0("UTR3", " (", gene_id, ")")
  } else if (type == "stop_codon") {
    anno <- paste0("Stop Codon", " (", gene_id, ")")
  }  else if (type == "start_codon") {
    anno <- paste0("Start Codon", " (", gene_id, ")")
  } else {
    anno <- paste0(type, " (", gene_id, ")")
  }
  dtt <- data.table::data.table(query_index = query_index, annotation = anno)
  dtt <- dtt[, toString(annotation), by = query_index]
  res <- list(query_index = query_index, annotation = anno, gene = gene_id)
  return(res)
}

.update_genomic_annotation <- function(peaks = NULL, genomic_region = NULL, type = NULL, anno = NULL, same_strand = FALSE) {
  hits <- .get_genomic_annotation_internal(peaks = peaks, genomic_region = genomic_region, type = type, same_strand = same_strand)
  if (length(hits) > 1) {
    hit_index <- hits$query_index
    anno[["annotation"]][hit_index] <- hits$annotation
    anno[["detail_genomic_annotation"]][hit_index, type] <- TRUE
  }
  return(anno)
}

.get_genomic_annotation <- function(peaks = NULL,  txdb_features = NULL, genomic_annotation_priority = NULL, same_strand = FALSE) {
  annotation <- rep(NA, length(peaks))
  flag <- rep(FALSE, length(peaks))
  detail_genomic_annotation <- data.frame(
    genic = flag,
    Intergenic = flag,
    Promoter = flag,
    five_utr = flag,
    three_utr = flag,
    CDS = flag,
    Intron = flag,
    stop_codon = flag)

  anno <- list(annotation = annotation, detail_genomic_annotation = detail_genomic_annotation)
  anno1 <- list(annotation = annotation, detail_genomic_annotation = detail_genomic_annotation)
  anno2 <- list(annotation = annotation, detail_genomic_annotation = detail_genomic_annotation)
  anno3 <- list(annotation = annotation, detail_genomic_annotation = detail_genomic_annotation)
  anno0 <- list(annotation = annotation, detail_genomic_annotation = detail_genomic_annotation)
  anno6 <- list(annotation = annotation, detail_genomic_annotation = detail_genomic_annotation)
  anno3 <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$exons, type = "CDS", anno = anno3, same_strand = same_strand)

  genomic_annotation_priority <- S4Vectors::rev(genomic_annotation_priority)
  for (AP in genomic_annotation_priority) {
    cat(paste0(AP,"\n"))
    if (AP == "Intergenic") {
      annotation[is.na(annotation)] <- "Intergenic"
      anno[["annotation"]] <- annotation
    } else if (AP == "Intron") {
      anno <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$introns, type = "Intron", anno = anno, same_strand = same_strand)
    } else if (AP == "CDS") {
      anno <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$cds, type = "CDS", anno = anno, same_strand = same_strand)
    } else if (AP == "3UTR") {
      anno <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$three_utrs, type = "three_utr", anno = anno, same_strand = same_strand)
    } else if (AP == "5UTR") {
      anno <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$five_utrs, type = "five_utr", anno = anno, same_strand = same_strand)
    } else if (AP == "stop_codon") {
      anno1 <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$stop_codon, type = "stop_codon", anno = anno1, same_strand = same_strand)
    } else if (AP == "Promoter") {
      anno <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$promoters, type = "Promoter", anno = anno, same_strand = same_strand)
    } else if (AP == "start_codon") {
      anno6 <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$start_codon, type = "start_codon", anno = anno6, same_strand = same_strand)
    } else if (AP == "Transcript") {
      anno0 <- .update_genomic_annotation(peaks = peaks, genomic_region = txdb_features$transcripts, type = "Transcript", anno = anno0, same_strand = same_strand)
    }
    annotation <- anno[["annotation"]]
    detail_genomic_annotation <- anno[["detail_genomic_annotation"]]
  }
  genic_index <- which(apply(detail_genomic_annotation[, c("CDS", "Intron","five_utr","three_utr")], 1, any))
  detail_genomic_annotation[-genic_index, "Intergenic"] <- TRUE
  detail_genomic_annotation[genic_index, "genic"] <- TRUE
  return(list(annotation = annotation,
              annotation_exon = anno3[["annotation"]],
              annotation_stop_codon = anno1[["annotation"]],
              annotation_start_codon = anno6[["annotation"]],
              annotation_transcript = anno0[["annotation"]],
              detail_genomic_annotation = detail_genomic_annotation))
}

.get_genomic_anno_stat <- function(peak_anno = NULL) {
  anno <- peak_anno$annotation
  anno_stop_codon  <- peak_anno$annotation_stop_codon
  anno_start_codon <- peak_anno$annotation_start_codon
  anno_transcript <- peak_anno$annotation_transcript
  e1 <- TRUE
  i1 <- TRUE
  ids <- TRUE
  if (is.null(e1) || !e1) {
    e1lab <- "1st CDS"
    anno[grep("CDS 1 of", anno)] <- e1lab
    exonlab <- "Other CDS"
  } else {
    e1lab <- NULL
    exonlab <- "CDS"
  }
  if (is.null(i1) || !i1) {
    i1lab <- "1st Intron"
    anno[grep("intron 1 of", anno)] <- i1lab
    intronlab <- "Other Intron"
  } else {
    i1lab <- NULL
    intronlab <- "Intron"
  }
  anno[grep("CDS \\(", anno)] <- exonlab
  anno[grep("Intron \\(", anno)] <- intronlab
  anno[grep("Promoter \\(", anno)] <- "Promoter"
  anno[grep("UTR5 \\(", anno)] <- "UTR5"
  anno[grep("UTR3 \\(", anno)] <- "UTR3"
  lvs <- c(
    "Promoter",
    "UTR5",
    "Start Codon",
    e1lab,
    exonlab,
    "Stop Codon",
    "UTR3",
    intronlab,
    "Intergenic"
  )
  anno_stop_codon <- anno_stop_codon[grep("Stop Codon \\(", anno_stop_codon)]
  anno_stop_codon[grep("Stop Codon \\(", anno_stop_codon)] <- "Stop Codon"
  anno_start_codon <- anno_start_codon[grep("^Start Codon", anno_start_codon)]
  anno_start_codon[grep("^Start Codon", anno_start_codon)] <- "Start Codon"
  anno_transcript <- anno_transcript[grep("^Transcript", anno_transcript)]
  anno_transcript[grep("^Transcript", anno_transcript)] <- "Transcript"
  anno_table <- table(anno)
  anno_stop_codon.table <- table(anno_stop_codon)
  anno_start_codon.table <- table(anno_start_codon)
  anno_transcript.table <- table(anno_transcript)
  anno_table <- c(anno_table, anno_start_codon.table, anno_stop_codon.table)
  anno_ratio <- anno_table
  anno_df <- data.table::as.data.table(anno_ratio, keep.rownames = TRUE)
  data.table::setnames(anno_df, c("Feature", "Frequency"))
  anno_df$Feature <- factor(anno_df$Feature, levels = lvs[lvs %in% anno_df$Feature])
  anno_df <- data.table::setDT(anno_df[order(anno_df$Feature), ])
  return(anno_df)
}

.get_peak_exon_stat <- function(peak_anno = NULL, tx_length_file = NULL) {
  anno <- peak_anno$annotation_transcript
  anno = anno[grep("^Transcript ", anno)]
  exon_rank = data.table::setDT(rtracklayer::as.data.frame(stringr::str_match(anno,"Transcript \\((\\S+)\\)")))
  exon_rank <- data.table::setnames(exon_rank,c("V2"), c("tx_name"))
  exon_rank <- IRanges::unique(exon_rank, by="tx_name")
  exon_rank = SummarizedExperiment::merge(exon_rank, tx_length_file, by="tx_name", all.x = TRUE)
  exon_rank = exon_rank[!is.na(tx_name)]
  return(exon_rank)
}

.get_peak_exon_size <- function(peak_anno = NULL, genomic_region = NULL, type = "CDS") {
  anno <- peak_anno$annotation_exon
  anno = anno[grep("^CDS ", anno)]
  exon_rank = data.table::setDT(rtracklayer::as.data.frame(stringr::str_match(anno, "CDS \\((\\S+),.*CDS (\\d+) of (\\d+), (\\d+)bp\\)")))
  exon_rank[, V3 := as.numeric(V3)]
  exon_rank[, V4 := as.numeric(V4)]
  exon_rank[, exonLoc := "Middle"]
  exon_rank[V3 == 1, exonLoc := "First"]
  exon_rank[V3==V4,exonLoc := "Last"]
  exon_rank[, exonID := paste0(V2, ":", V3)]
  data.table::setnames(exon_rank, c("V2", "V3"), c("tx_name", "rank"))
  if(class(genomic_region) == "GRangesList" | class(genomic_region) == "CompressedGRangesList"){
    GRegion <- BiocGenerics::unlist(genomic_region)
    GRegion$tx_name <- names(GRegion)
    GRegion$length = rtracklayer::width(GRegion)
  } else {
    GRegion <- genomic_region
    GRegion$length = rtracklayer::width(GRegion)
  }
  if (type == "Intron" || type =="CDS") {
    gr2 <- GRegion[!SummarizedExperiment::duplicated(GRegion$tx_name)]
    temp = data.table::data.table(name = GRegion$tx_name)
    temp = temp[, .N, by = name]
    strd <- as.character(rtracklayer::strand(gr2))
    GRegion_length <- temp[, N]
    names(GRegion_length) <- temp[, name]
    GRegion$rank <- lapply(seq_along(strd), function(i) {
      rank <- seq(1, GRegion_length[i])
      if (strd[i] == '-')
        rank <- rev(rank)
      return(rank)
    }) %>% unlist
  }

  exon_info = data.table::setDT(rtracklayer::as.data.frame(GenomicRanges::mcols(GRegion)))
  exon_info[, exonID := paste0(tx_name, ":", rank)]
  exon_info = SummarizedExperiment::merge(exon_info, exon_rank, by="exonID", all.x=T)
  exon_info = exon_info[!is.na(exonLoc)]
  exon_info[is.na(exonLoc), exonLoc := "Other"]
  return(exon_info)
}

.get_gene_matrix <- function(peak_anno = NULL, peak_obj = NULL, gff_file = NULL) {
  anno <- peak_anno$annotation
  anno_stop_codon <- peak_anno$annotation_stop_codon
  anno_transcript <- peak_anno$annotation_transcript
  peak_opr <- data.table::copy(peak_obj)
  peak_opr[, Index := .I]
  peak_opr[, peak := paste0(seqnames, ":", start, "-", end, "|", peak_name)]

  exon_rank <- data.table::setDT(rtracklayer::as.data.frame(stringr::str_match(anno, "(.*) \\((.+)\\)")))
  data.table::setnames(exon_rank, c("V2","V3"), c("feature","transcript_id"))
  exon_rank[,transcript_id := IRanges::gsub(",.*$","",transcript_id)]
  exon_rank2 <- data.table::copy(exon_rank)
  exon_rank2[, Index := .I]
  exon_rank2 <- SummarizedExperiment::merge(exon_rank2,peak_opr[ ,c("peak","Index")], by = "Index", all.x = TRUE)
  exon_rank2 <- exon_rank2[!is.na(feature)]
  exon_rank <- exon_rank[!is.na(feature)]

  stop_codon_rank <- data.table::setDT(rtracklayer::as.data.frame(stringr::str_match(anno_stop_codon, "(.*) \\((.+)\\)")))
  data.table::setnames(stop_codon_rank, c("V2", "V3"), c("feature", "transcript_id"))
  stop_codon_rank <- stop_codon_rank[!is.na(feature)]
  transcript_rank <- data.table::setDT(as.data.frame(str_match(anno_transcript, "(.*) \\((.+)\\)")))
  setnames(transcript_rank, c("V2", "V3"), c("feature", "transcript_id"))
  transcript_rank <- transcript_rank[!is.na(feature)]
  exon_rank <- data.table::rbindlist(list(exon_rank, stop_codon_rank, transcript_rank))
  exon_rank <- exon_rank[ , .N, by = c('transcript_id', 'feature')]
  peak_detail <- exon_rank2[, list(peak_detail = paste(peak, collapse = ",")), by = transcript_id]
  feature_stat <- data.table::setDT(data.table::dcast(exon_rank, transcript_id ~ feature, value.var = "N"))
  feature_stat[is.na(feature_stat)] <- 0
  obj_names <- c("UTR5", "CDS", "UTR3", "Intron")
  obj_names <- obj_names[obj_names %in% names(feature_stat)]

  feature_stat[,Transcript := MatrixGenerics::rowSums(.SD), .SDcols = obj_names]
  cct <- data.table::setDT(rtracklayer::as.data.frame(gff_file))
  cct <- cct[type == "transcript"]
  results <- merge(feature_stat, cct[ , c("gene_id", "transcript_id", "gene_type", "gene_name")], by="transcript_id", all.x = TRUE)
  results <- merge(results, peak_detail, by="transcript_id", all.x = TRUE)
  order <- c("gene_id","transcript_id", "gene_name", "gene_type", "Promoter", "UTR5", "CDS", "UTR3", "Stop Codon", "Intron", "Transcript","peak_detail")
  for(i in order) {
    if(! i %in% names(results)) {
      results[, S4Vectors::eval(i) := 0]
    }
  }
  data.table::setcolorder(results, order)
  return(results)
}

.batch_genomic_annotation <- function(index = NULL,
                                      peak_gr = NULL,
                                      txdb_features = NULL,
                                      genomic_annotation_priority = NULL,
                                      tx_length_file = NULL,
                                      gff_file = NULL,
                                      peak_dt = NULL,
                                      group_name = NULL) {
  peak_anno <- .get_genomic_annotation(peaks = peak_gr[[index]], txdb_features = txdb_features, genomic_annotation_priority = genomic_annotation_priority)
  rp <- .get_genomic_anno_stat(peak_anno = peak_anno)
  rp1 <- .get_peak_exon_stat(peak_anno = peak_anno, tx_length_file = tx_length_file)
  rp2 <- .get_peak_exon_size(peak_anno = peak_anno, genomic_region = txdb_features$exons)
  rp3 <- .get_gene_matrix(peak_anno = peak_anno, peak_obj = peak_dt[[index]], gff_file = gff_file)
  rp6 <- data.table::data.table(index = 1:length(peak_anno$annotation_transcript), gene = peak_anno$annotation_transcript, anno = peak_anno$annotation)
  temp_dt <- data.table::copy(peak_dt[[index]])
  temp_dt[, peak := paste0(seqnames, ":", start, "-", end,"|", peak_name)]
  temp_dt[, index := .I]
  rp6 <- rp6[temp_dt[,list(seqnames, start, end, width, index, peak, group_name)], on = "index"]
  rp6[!is.na(gene), gene := anno]
  rp6 <- rp6[!is.na(gene)]

  rp6[, transcript_id := IRanges::gsub(".*\\((\\S+).*\\)", "\\1", gene, perl = TRUE)]
  rp6[, transcript_id := IRanges::gsub(",", "", transcript_id)]
  cct <- data.table::setDT(rtracklayer::as.data.frame(gff_file))
  cct <- cct[type == "transcript"]
  rp6 <- SummarizedExperiment::merge(rp6, cct[ , c("transcript_id", "gene_type", "gene_name")], by = "transcript_id", all.x = TRUE)
  rp6[,location := IRanges::gsub("(.*) \\(\\S+.*","\\1", anno, perl = TRUE)]
  rp6[location == "Intergenic", location := "Transcript"]
  rp6[, transcriptID := transcript_id]
  rp[, Sample := names(peak_gr[index])]
  rp1[, Sample := names(peak_gr[index])]
  rp2[, Sample := names(peak_gr[index])]
  rp3[, Sample := names(peak_gr[index])]
  gn <- names(peak_gr[index])
  rp6[, Sample := gn]
  rp6[, c("transcript_id", "index", "gene") := NULL]
  res <- list(rp, rp2, rp3, rp6, rp1)
  return(res)
}
