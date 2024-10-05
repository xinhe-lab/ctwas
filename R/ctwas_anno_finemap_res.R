
#' Get original molecular IDs
#'
#' @param finemap_res a data frame of cTWAS finemapping results.
#'
#' @return a vector of original molecular IDs
#'
#' @export
get_molecular_ids <- function(finemap_res){
  molecular_ids <- finemap_res$id
  gidx <- which(finemap_res$group!="SNP")
  gids <- finemap_res$id[gidx]
  molecular_ids[gidx] <- sapply(strsplit(gids, split = "[|]"), "[[", 1)
  return(molecular_ids)
}

#' convert z-scores to pvalues
#'
#' @param z a vector of z-scores
#'
#' @param neg_log10_p If TRUE, returns -log10(p-values). Otherwise, returns p-values.
#'
#' @return a vector of p-values, or -log10(p-values) if neg_log10_p is TRUE.
#'
#' @export
z2p <- function(z, neg_log10_p = FALSE){
  if (neg_log10_p) {
    nlog10p <- -(log(2) + pnorm(abs(z), lower.tail=FALSE, log.p=TRUE))/log(10)
    return(nlog10p)
  } else {
    p <- 2*pnorm(abs(z), lower.tail=FALSE)
    return(p)
  }
}

#' @title Map finemapping result of molecular traits to genes.
#'
#' @param finemap_res a data frame of cTWAS finemapping results.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param map_by column name to be mapped by (default: "molecular_id").
#'
#' @param add_gene_annot If TRUE, add annotations
#'
#' @param add_position If TRUE, add positions
#'
#' @param use_gene_pos Use mid (midpoint), start, or end positions as gene positions.
#'
#' @param drop_unmapped If TRUE, remove unmapped genes.
#'
#' @return a data frame of annotated finemapping result.
#'
#' @importFrom logging loginfo
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate select group_by ungroup n
#' @importFrom rlang .data
#'
#' @export
#'
anno_finemap_res <- function(finemap_res,
                             snp_map,
                             mapping_table,
                             map_by = "molecular_id",
                             add_gene_annot = TRUE,
                             add_position = TRUE,
                             use_gene_pos = c("mid", "start", "end"),
                             drop_unmapped = TRUE){

  loginfo("Annotating fine-mapping result ...")

  use_gene_pos <- match.arg(use_gene_pos)

  # map molecular traits to genes
  # there could be n-to-1 and 1-to-n mapping
  if (add_gene_annot) {
    if (!is.null(finemap_res$gene_name)){
      loginfo("'gene_name' is already available. Skip annotating genes.")
    } else {
      if (!map_by %in% colnames(finemap_res)) {
        stop(paste("column", map_by, "cannot be found in finemap_res!"))
      }

      if (!map_by %in% colnames(mapping_table)) {
        stop(paste("column", map_by, "cannot be found in mapping_table!"))
      }

      # gene results
      finemap_gene_res <- finemap_res[finemap_res$group!="SNP",]

      # extract gene ids
      if (is.null(finemap_gene_res$molecular_id)) {
        finemap_gene_res$molecular_id <- sapply(strsplit(finemap_gene_res$id, split = "[|]"), "[[", 1)
      }

      # Map molecular traits to genes by joining finemap_gene_res with mapping_table
      loginfo("Map molecular traits to genes")
      finemap_gene_res <- finemap_gene_res %>%
        left_join(mapping_table, by = map_by, multiple = "all")

      unmapped_finemap_res <- NULL
      if (drop_unmapped) {
        unmapped_idx <- which(is.na(finemap_gene_res$gene_name))
        if ( length(unmapped_idx) > 0){
          unmapped_molecular_ids <- unique(finemap_gene_res$molecular_id[unmapped_idx])
          loginfo("Drop %d unmapped molecular traits", length(unmapped_molecular_ids))
          finemap_gene_res <- finemap_gene_res[-unmapped_idx, , drop=FALSE]
        }
      }

      # split PIPs for molecular traits (e.g. introns) mapped to multiple genes
      if (any(duplicated(finemap_gene_res$id))) {
        loginfo("Split PIPs for molecular traits mapped to multiple genes")
        finemap_gene_res <- finemap_gene_res %>%
          group_by(.data$id) %>%
          mutate(susie_pip = ifelse(n() > 1, .data$susie_pip / n(), .data$susie_pip)) %>%
          ungroup()
        finemap_gene_res <- as.data.frame(finemap_gene_res)
      }

      # SNP results
      finemap_snp_res <- finemap_res[finemap_res$group=="SNP",]
      added_colnames <- setdiff(colnames(finemap_gene_res), colnames(finemap_snp_res))
      finemap_snp_res[, added_colnames] <- NA

      finemap_res <- rbind(finemap_gene_res, finemap_snp_res)
    }

  }

  # add positions
  if (add_position) {
    finemap_res <- add_pos(finemap_res,
                           snp_map = snp_map,
                           mapping_table = mapping_table,
                           map_by = map_by,
                           use_gene_pos = use_gene_pos)
  }


  return(finemap_res)
}

#' @title Map molecular traits to genes in susie alpha result.
#'
#' @param susie_alpha_res a data frame of susie alpha result.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "molecular_id", "gene_name".
#'
#' @param map_by column name to be mapped by (default: "molecular_id").
#'
#' @param drop_unmapped If TRUE, remove unmapped genes.
#'
#' @return a data frame of annotated finemapping result.
#'
#' @importFrom logging loginfo
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate select group_by ungroup n
#' @importFrom rlang .data
#'
#' @export
#'
anno_susie_alpha_res <- function(susie_alpha_res,
                                 mapping_table,
                                 map_by = "molecular_id",
                                 drop_unmapped = TRUE){

  loginfo("Annotating susie alpha result ...")

  # map molecular traits to genes
  # there could be n-to-1 and 1-to-n mapping
  if (!is.null(susie_alpha_res$gene_name)){
    loginfo("'gene_name' is already available. Skip annotating genes.")
  } else {
    if (!map_by %in% colnames(susie_alpha_res)) {
      stop(paste("column", map_by, "cannot be found in susie_alpha_res!"))
    }

    if (!map_by %in% colnames(mapping_table)) {
      stop(paste("column", map_by, "cannot be found in mapping_table!"))
    }

    # gene results
    susie_alpha_res <- susie_alpha_res[susie_alpha_res$group!="SNP",]

    # extract molecular ids
    if (is.null(susie_alpha_res$molecular_id)) {
      susie_alpha_res$molecular_id <- sapply(strsplit(susie_alpha_res$id, split = "[|]"), "[[", 1)
    }

    # Map molecular traits to genes by joining finemap_gene_res with mapping_table
    loginfo("Map molecular traits to genes")
    susie_alpha_res <- susie_alpha_res %>%
      left_join(mapping_table, by = map_by, multiple = "all")

    if (drop_unmapped) {
      unmapped_idx <- which(is.na(susie_alpha_res$gene_name))
      if ( length(unmapped_idx) > 0){
        unmapped_molecular_ids <- unique(susie_alpha_res$molecular_id[unmapped_idx])
        loginfo("Drop %d unmapped molecular traits", length(unmapped_molecular_ids))
        susie_alpha_res <- susie_alpha_res[-unmapped_idx, , drop=FALSE]
      }
    }

    # add a column with both ID and susie set information
    susie_alpha_res$tmp_id <- paste0(susie_alpha_res$id, ".", susie_alpha_res$susie_set)

    # split PIPs for molecular traits (e.g. introns) mapped to multiple genes
    if (any(duplicated(susie_alpha_res$tmp_id))) {
      loginfo("Split PIPs for molecular traits mapped to multiple genes")
      susie_alpha_res <- susie_alpha_res %>%
        group_by(.data$tmp_id) %>%
        mutate(susie_alpha = ifelse(n() > 1, .data$susie_alpha / n(), .data$susie_alpha)) %>%
        ungroup()
    }

    susie_alpha_res <- susie_alpha_res[, colnames(susie_alpha_res)!="tmp_id"]
    susie_alpha_res <- as.data.frame(susie_alpha_res)

  }

  return(susie_alpha_res)
}



#' @title Add SNP and gene positions to cTWAS finemapping result.
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param mapping_table a data frame of gene annotations, with required columns:
#' "gene_name", "chrom", "pos" (or "start" and "end").
#'
#' @param map_by column name to be mapped by (default: "gene_name").
#'
#' @param use_gene_pos Use mid (midpoint), start, or end positions as gene positions.
#'
#' @return a data frame with chromosomes and positions added to fine-mapping result
#'
#' @importFrom stats na.omit
#' @importFrom logging loginfo
#' @importFrom data.table rbindlist
#' @importFrom readr parse_number
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate select group_by ungroup n
#' @importFrom rlang .data
#'
#' @export
#'
add_pos <- function(finemap_res,
                    snp_map,
                    mapping_table,
                    map_by = "molecular_id",
                    use_gene_pos = c("mid", "start", "end")){

  # add gene positions
  loginfo("Add gene positions")
  finemap_gene_res <- finemap_res[finemap_res$group!="SNP",]
  gene_idx <- match(finemap_gene_res[,map_by], mapping_table[,map_by])
  tmp_mapping_table <- mapping_table[gene_idx, ]
  tmp_mapping_table$chrom <- parse_number(as.character(tmp_mapping_table$chrom))
  tmp_mapping_table$start <- as.numeric(tmp_mapping_table$start)
  tmp_mapping_table$end <- as.numeric(tmp_mapping_table$end)
  if (use_gene_pos == "mid"){
    finemap_gene_res$pos <- round((tmp_mapping_table$start+tmp_mapping_table$end)/2)
  } else if (use_gene_pos == "start") {
    finemap_gene_res$pos <- tmp_mapping_table$start
  } else if (use_gene_pos == "end") {
    finemap_gene_res$pos <- tmp_mapping_table$end
  }

  # add SNP positions
  loginfo("Add SNP positions")
  finemap_snp_res <- finemap_res[finemap_res$group=="SNP",]
  snp_info <- as.data.frame(rbindlist(snp_map))
  snp_idx <- match(finemap_snp_res$id, snp_info$id)
  finemap_snp_res$chrom <- snp_info$chrom[snp_idx]
  finemap_snp_res$chrom <- parse_number(as.character(finemap_snp_res$chrom))
  finemap_snp_res$pos <- as.numeric(snp_info$pos[snp_idx])

  new_colnames <- unique(c(colnames(finemap_res), "chrom", "pos"))

  return(rbind(finemap_gene_res[,new_colnames],
               finemap_snp_res[,new_colnames]))
}

#' @title Get gene annotation table from Ensembl database
#' @param ens_db Ensembl gene annotation database
#' @param gene_ids Ensembl gene IDs
#'
#' @importFrom stats na.omit
#' @importFrom ensembldb genes
#' @importFrom AnnotationFilter GeneIdFilter
#' @importFrom logging loginfo
#'
#' @export
get_gene_annot_from_ens_db <- function(ens_db, gene_ids) {

  gene_ids <- unique(na.omit(gene_ids))
  if (any(grep("[.]", gene_ids))) {
    gene_ids_trimmed <- sapply(strsplit(gene_ids, split = "[.]"), "[[", 1)
  }
  gene_annot_gr <- genes(ens_db, filter = GeneIdFilter(gene_ids_trimmed))
  annot_idx <- match(gene_ids_trimmed, gene_annot_gr$gene_id)
  if (anyNA(annot_idx)) {
    loginfo("Remove %d gene_ids not found in ens_db.", length(which(is.na(annot_idx))))
    gene_ids <- gene_ids[!is.na(annot_idx)]
    gene_ids_trimmed <- gene_ids_trimmed[!is.na(annot_idx)]
    annot_idx <- annot_idx[!is.na(annot_idx)]
  }
  gene_annot <- as.data.frame(gene_annot_gr[annot_idx,])
  gene_annot$gene_id <- gene_ids
  colnames(gene_annot)[colnames(gene_annot) == 'seqnames'] <- 'chrom'
  colnames(gene_annot)[colnames(gene_annot) == 'gene_biotype'] <- 'gene_type'
  gene_annot$start <- as.numeric(gene_annot$start)
  gene_annot$end <- as.numeric(gene_annot$end)
  gene_annot <- gene_annot[, c("gene_id", "gene_name", "gene_type", "chrom", "start", "end")]

  return(gene_annot)
}

