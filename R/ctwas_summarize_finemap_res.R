
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
          # unmapped_finemap_res <- finemap_gene_res[unmapped_idx, ]
          loginfo("Drop %d unmapped molecular IDs", length(unmapped_idx))
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

#' @title Combines gene PIPs by context, type or group.
#'
#' @param finemap_res a data frame of annotated cTWAS finemapping result
#'
#' @param group_by column name to group genes by.
#'
#' @param by option to combine PIPs by: "context" (default), "type", or "group".
#'
#' @param method method to combine PIPs of molecular traits targeting the same gene.
#' options:
#' "combine_cs" (default): first sum PIPs of molecular traits for the same gene in the same CS,
#' then apply the multiplication formula across CS.
#' "sum" sum over PIPs of all molecular traits for the same gene;
#' "combine_all" use the multiplication formula for all molecular traits for the same gene.
#'
#' @param filter_cs If TRUE, limits gene results to credible sets.
#'
#' @param missing_value set missing value as (default: NA)
#'
#' @return a data frame of combined gene PIPs for each context, type or group
#'
#' @importFrom magrittr %>%
#' @importFrom stats aggregate
#' @importFrom dplyr left_join
#'
#' @export
combine_gene_pips <- function(finemap_res,
                              group_by = "gene_name",
                              by = c("context", "type", "group"),
                              method = c("combine_cs", "combine_all", "sum"),
                              filter_cs = TRUE,
                              missing_value = NA){

  by <- match.arg(by)
  method <- match.arg(method)

  # Check to see if gene_name and gene_type are already in finemap_res
  if (!(group_by %in% colnames(finemap_res))){
    stop(paste("Cannot find the column", group_by, "in finemap_res!"))
  }

  # cannot filter CS or combine CS for finemapping results from "no-LD" version
  if (is.null(finemap_res$cs_index)){
    filter_cs <- FALSE
    method <- "sum"
  }

  # work with gene results below
  finemap_gene_res <- finemap_res[finemap_res$group!="SNP",]

  # limit genes to credible sets
  if (filter_cs) {
    loginfo("Limit gene results to credible sets")
    finemap_gene_res <- finemap_gene_res[finemap_gene_res$cs_index!=0,]
  }

  # combine PIPs across all contexts or types
  combined_gene_pips <- compute_combined_pips(finemap_gene_res,
                                              group_by = group_by,
                                              method = method,
                                              filter_cs = filter_cs)
  colnames(combined_gene_pips) <- c(group_by, "combined_pip")

  if (by == "context"){
    # combine PIPs for each context
    contexts <- unique(finemap_gene_res$context)
    for (context in contexts){
      tmp_finemap_gene_res <- finemap_gene_res[finemap_gene_res$context==context,,drop=FALSE]
      tmp_combined_gene_pips <- compute_combined_pips(tmp_finemap_gene_res,
                                                      group_by = group_by,
                                                      method = method,
                                                      filter_cs = filter_cs)
      colnames(tmp_combined_gene_pips) <- c(group_by, paste0(context, "_pip"))
      combined_gene_pips <- combined_gene_pips %>%
        left_join(tmp_combined_gene_pips, by = group_by)
    }
  } else if (by == "type") {
    # combine PIPs for each type
    types <- unique(finemap_gene_res$type)
    for (type in types){
      tmp_finemap_gene_res <- finemap_gene_res[finemap_gene_res$type==type,,drop=FALSE]
      tmp_combined_gene_pips <- compute_combined_pips(tmp_finemap_gene_res,
                                                      group_by = group_by,
                                                      method = method,
                                                      filter_cs = filter_cs)
      colnames(tmp_combined_gene_pips) <- c(group_by, paste0(type, "_pip"))
      combined_gene_pips <- combined_gene_pips %>%
        left_join(tmp_combined_gene_pips, by = group_by)
    }
  } else if (by == "group") {
    # combine PIPs for each group
    groups <- unique(finemap_gene_res$group)
    for (group in groups){
      tmp_finemap_gene_res <- finemap_gene_res[finemap_gene_res$group==group,,drop=FALSE]
      tmp_combined_gene_pips <- compute_combined_pips(tmp_finemap_gene_res,
                                                      group_by = group_by,
                                                      method = method,
                                                      filter_cs = filter_cs)
      colnames(tmp_combined_gene_pips) <- c(group_by, paste0(group, "_pip"))
      combined_gene_pips <- combined_gene_pips %>%
        left_join(tmp_combined_gene_pips, by = group_by)
    }
  }

  if (!is.na(missing_value)) {
    combined_gene_pips[is.na(combined_gene_pips)] <- missing_value
  }

  # order by combined PIP
  combined_gene_pips <- combined_gene_pips[order(-combined_gene_pips$combined_pip),]
  rownames(combined_gene_pips) <- NULL

  new_colnames <- c(setdiff(colnames(combined_gene_pips), "combined_pip"), "combined_pip")
  combined_gene_pips <- combined_gene_pips[, new_colnames]

  return(combined_gene_pips)
}

# Compute combined gene PIP using the multiplication formula
# combined gene PIP = 1 - \prod_k (1 - PIP_k).
# PIP_k is the PIP of the k-th molecular trait of a gene.
combine_pips <- function(pips){
  return(1 - prod(1 - pips))
}

# Computes combined gene PIPs.
# "combine_cs" (default): first sum PIPs of molecular traits for the same gene in the same CS,
# then apply the multiplication formula across CS.
# "sum" sum over PIPs of all molecular traits for the same gene;
# "combine_all" use the multiplication formula for all molecular traits for the same gene.
compute_combined_pips <- function(finemap_gene_res,
                                  group_by = "gene_name",
                                  method = c("combine_cs", "combine_all", "sum"),
                                  filter_cs = TRUE){

  method <- match.arg(method)

  if (filter_cs) {
    finemap_gene_res <- finemap_gene_res[finemap_gene_res$cs_index!=0,]
  }

  if (method == "sum") {
    combined_gene_pips <- aggregate(data.frame(susie_pip = finemap_gene_res$susie_pip),
                                    by = list(gene = finemap_gene_res[,group_by]),
                                    FUN = sum)
    colnames(combined_gene_pips) <- c(group_by, "combined_pip")
  } else if (method == "combine_all") {
    combined_gene_pips <- aggregate(data.frame(susie_pip = finemap_gene_res$susie_pip),
                                    by = list(gene = finemap_gene_res[,group_by]),
                                    FUN = combine_pips)
    colnames(combined_gene_pips) <- c(group_by, "combined_pip")
  } else if (method == "combine_cs") {
    finemap_gene_res$cs_id <- paste0(finemap_gene_res$region_id, ".", finemap_gene_res$cs_index)
    # for each gene, first sum PIPs within the same credible sets
    summed_gene_pips <- aggregate(data.frame(susie_pip = finemap_gene_res$susie_pip),
                                  by = list(gene = finemap_gene_res[,group_by],
                                            cs_id = finemap_gene_res[,"cs_id"]),
                                  FUN = sum)
    # then combine PIPs using the multiplication formula across credible sets
    combined_gene_pips <- aggregate(data.frame(susie_pip = summed_gene_pips$susie_pip),
                                    by = list(gene = summed_gene_pips$gene),
                                    FUN = combine_pips)
    colnames(combined_gene_pips) <- c(group_by, "combined_pip")
  }

  return(combined_gene_pips)
}

