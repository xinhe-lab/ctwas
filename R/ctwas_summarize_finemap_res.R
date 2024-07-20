
#' @title Annotates cTWAS finemapping result with gene annotations
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param snp_map a list of SNP-to-region map for the reference.
#'
#' @param gene_annot a data frame of gene annotations, with columns:
#' "chrom", "start", "end", "gene_id", "gene_name", "gene_type.
#'
#' @param use_gene_pos use mid (midpoint), start or end positions to
#' represent gene positions.
#'
#' @param drop_unannotated_genes If TRUE, remove unannotated genes.
#'
#' @param filter_protein_coding_genes If TRUE, keep protein coding genes only.
#'
#' @param filter_cs If TRUE, limits results in credible sets.
#'
#' @return a data frame of cTWAS finemapping result including gene
#' names, types and positions.
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
anno_finemap_res <- function(finemap_res,
                             snp_map,
                             gene_annot = NULL,
                             use_gene_pos = c("mid", "start", "end"),
                             drop_unannotated_genes = TRUE,
                             filter_protein_coding_genes = FALSE,
                             filter_cs = FALSE){

  loginfo("Annotating ctwas fine-mapping result ...")

  use_gene_pos <- match.arg(use_gene_pos)

  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))

  # Check required columns for gene_annot
  annot_cols <- c("gene_id", "gene_name", "gene_type", "start", "end")
  if (!all(annot_cols %in% colnames(gene_annot))){
    stop("gene_annot needs to contain the following columns: ",
         paste(annot_cols, collapse = " "))
  }

  # limit to protein coding genes
  if (filter_protein_coding_genes) {
    loginfo("Limit to protein coding genes")
    gene_annot <- gene_annot[gene_annot$gene_type=="protein_coding",]
  }

  # limit to credible sets
  if (filter_cs) {
    loginfo("Limit to results in credible sets")
    finemap_res <- finemap_res[finemap_res$cs_index!=0,]
  }

  # extract gene ids
  finemap_gene_res <- finemap_res[finemap_res$type!="SNP",]

  if (is.null(finemap_gene_res$gene_id)) {
    finemap_gene_res$gene_id <- sapply(strsplit(finemap_gene_res$id, split = "[|]"), "[[", 1)
  }

  # add gene_name and gene_type
  if (!all(c("gene_name", "gene_type") %in% colnames(finemap_gene_res))) {
    # avoid duplicating columns when doing left_join()
    finemap_gene_res[, c("gene_name", "gene_type", "start", "end")] <- NULL
    if (!is.null(gene_annot$chrom)) {
      finemap_gene_res$chrom <- NULL
    }

    # add gene annotations
    loginfo("Add gene annotations")
    finemap_gene_res <- finemap_gene_res %>%
      left_join(gene_annot, by = "gene_id", multiple = "all") %>%
      mutate(start = as.numeric(.data$start), end = as.numeric(.data$end)) %>%
      mutate(chrom = parse_number(as.character(.data$chrom)))

    if (drop_unannotated_genes) {
      loginfo("Remove unannotated genes")
      finemap_gene_res <- na.omit(finemap_gene_res)
    }

    # split PIPs for molecular traits (e.g. introns) mapped to multiple genes
    if (any(duplicated(finemap_gene_res$id))) {
      loginfo("Split PIPs for traits mapped to multiple genes")
      finemap_gene_res <- finemap_gene_res %>%
        group_by(.data$id) %>%
        mutate(susie_pip = ifelse(n() > 1, .data$susie_pip / n(), .data$susie_pip)) %>%
        ungroup()
    }
  }

  # update gene position info
  loginfo("Use gene %s positions", use_gene_pos)
  if (use_gene_pos == "mid"){
    # use the midpoint as gene position
    finemap_gene_res$pos <- round((finemap_gene_res$start+finemap_gene_res$end)/2)
  } else if (use_gene_pos == "start") {
    finemap_gene_res$pos <- finemap_gene_res$start
  } else if (use_gene_pos == "end") {
    finemap_gene_res$pos <- finemap_gene_res$end
  }

  # add SNP positions
  loginfo("Add SNP positions")
  finemap_snp_res <- finemap_res[finemap_res$type=="SNP",]
  snp_idx <- match(finemap_snp_res$id, snp_info$id)
  finemap_snp_res$chrom <- snp_info$chrom[snp_idx]
  finemap_snp_res$chrom <- parse_number(as.character(finemap_snp_res$chrom))
  finemap_snp_res$pos <- snp_info$pos[snp_idx]
  finemap_snp_res$gene_id <- NA
  finemap_snp_res$gene_name <- NA
  finemap_snp_res$gene_type <- "SNP"

  new_colnames <- unique(c("chrom", "pos", "gene_id", "gene_name", "gene_type", colnames(finemap_res)))
  finemap_res <- rbind(finemap_gene_res[,new_colnames],
                       finemap_snp_res[,new_colnames])

  return(finemap_res)
}


#' @title Get gene annotation table from Ensembl database
#' @param ens_db Ensembl gene annotation database
#' @param gene_ids Ensembl gene IDs
#'
#' @importFrom stats na.omit
#' @importFrom ensembldb genes
#' @importFrom AnnotationFilter GeneIdFilter
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
    warning(sprintf("remove %d gene_ids not found in ens_db.",
                    length(which(is.na(annot_idx)))))
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
  gene_annot <- gene_annot[, c("chrom", "start", "end", "gene_id", "gene_name", "gene_type")]
  return(gene_annot)
}


#' @title Combines gene PIPs by context, type or group
#'
#' @param finemap_res a data frame of annotated cTWAS finemapping result
#'
#' @param by sum gene PIPs by "context", "type", or "group".
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only.
#'
#' @param filter_cs TRUE/FALSE. If TRUE, limits results in credible sets.
#'
#' @param replace_NA values to replace NAs with.
#'
#' @param digits digits to round combined PIPs
#'
#' @return a data frame of combined gene PIPs for each context, type or group
#'
#' @importFrom magrittr %>%
#' @importFrom stats aggregate
#' @importFrom dplyr left_join
#'
#' @export
combine_gene_pips <- function(finemap_res,
                              by = c("context", "type", "group"),
                              filter_protein_coding_genes = FALSE,
                              filter_cs = FALSE,
                              replace_NA = NA,
                              digits = 3){

  by <- match.arg(by)

  # Check to see if gene_name and gene_type are already in finemap_res
  annot_cols <- c("gene_name", "gene_type")
  if (!all(annot_cols %in% colnames(finemap_res))){
    stop("finemap_res needs to contain the following columns: ",
         paste(annot_cols, collapse = " "),
         "\nPlease first run anno_finemap_res() to annotate finemap_res")
  }

  finemap_gene_res <- finemap_res[finemap_res$type!="SNP",]

  # limit to protein coding genes
  if (filter_protein_coding_genes) {
    finemap_gene_res <- finemap_gene_res[finemap_gene_res$gene_type=="protein_coding",]
  }

  # filter credible sets
  if (filter_cs) {
    finemap_gene_res <- finemap_gene_res[finemap_gene_res$cs_index!=0,]
  }

  # sum PIPs
  combined_gene_pips <- aggregate(finemap_gene_res$susie_pip,
                                  by=list(finemap_gene_res$gene_name),
                                  FUN=sum)
  colnames(combined_gene_pips) <- c("gene_name", "combined_pip")

  if (by == "context"){
    # sum PIPs for each context
    contexts <- unique(finemap_gene_res$context)
    for (context in contexts){
      tmp_res <- finemap_gene_res[finemap_gene_res$context==context, c("gene_name", "susie_pip")]
      tmp_res <- aggregate(tmp_res$susie_pip, by=list(tmp_res$gene_name), FUN=sum)
      colnames(tmp_res) <- c("gene_name", paste0(context, "_pip"))
      combined_gene_pips <- combined_gene_pips %>% left_join(tmp_res, by = "gene_name")
    }
  } else if (by == "type") {
    # sum PIPs for each type
    types <- unique(finemap_gene_res$type)
    for (type in types){
      tmp_res <- finemap_gene_res[finemap_gene_res$type==type, c("gene_name", "susie_pip")]
      tmp_res <- aggregate(tmp_res$susie_pip, by=list(tmp_res$gene_name), FUN=sum)
      colnames(tmp_res) <- c("gene_name", paste0(type, "_pip"))
      combined_gene_pips <- combined_gene_pips %>% left_join(tmp_res, by = "gene_name")
    }
  } else if (by == "group") {
    # sum PIPs for each group
    groups <- unique(finemap_gene_res$group)
    for (group in groups){
      tmp_res <- finemap_gene_res[finemap_gene_res$group==group, c("gene_name", "susie_pip")]
      tmp_res <- aggregate(tmp_res$susie_pip, by=list(tmp_res$gene_name), FUN=sum)
      colnames(tmp_res) <- c("gene_name", paste0(group, "_pip"))
      combined_gene_pips <- combined_gene_pips %>% left_join(tmp_res, by = "gene_name")
    }
  }

  if (!is.na(replace_NA)) {
    combined_gene_pips[is.na(combined_gene_pips)] <- replace_NA
  }

  # order by combined PIP
  combined_gene_pips <- combined_gene_pips[order(-combined_gene_pips$combined_pip),]

  # round gene PIPs
  combined_gene_pips[,-1] <- round(combined_gene_pips[, -1], digits)

  return(combined_gene_pips)
}

