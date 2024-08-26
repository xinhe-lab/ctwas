
#' @title Map finemapping result of molecular traits to genes.
#'
#' @param finemap_res a data frame of cTWAS finemapping results.
#'
#' @param mapping_table a data frame of mapping between molecular traits and genes,
#' with required columns: "gene_id", "gene_name".
#'
#' @param map_by column name to be mapped by (default: "gene_id").
#'
#' @param drop_unmapped_genes If TRUE, remove unmapped genes.
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
                             mapping_table = NULL,
                             map_by = "gene_id",
                             drop_unmapped_genes = TRUE){

  loginfo("Mapping molecular traits to genes ...")

  # gene results
  finemap_gene_res <- finemap_res[finemap_res$group!="SNP",]

  # extract gene ids
  if (is.null(finemap_gene_res$gene_id)) {
    finemap_gene_res$gene_id <- sapply(strsplit(finemap_gene_res$id, split = "[|]"), "[[", 1)
  }

  # Check required columns for mapping_table
  required_cols <- c(map_by, "gene_name")
  if (!all(required_cols %in% colnames(mapping_table))){
    stop("mapping_table needs to contain the following columns: ",
         paste(required_cols, collapse = " "))
  }

  # map molecular traits to genes
  loginfo("Map molecular traits to genes")

  # there could be n-to-1 and 1-to-n mapping
  finemap_gene_res <- finemap_gene_res %>%
    left_join(mapping_table, by = map_by, multiple = "all")

  if (drop_unmapped_genes) {
    if ( any(is.na(finemap_gene_res$gene_name)) ){
      loginfo("Drop unmapped genes")
      finemap_gene_res <- finemap_gene_res[!is.na(finemap_gene_res$gene_name), , drop=FALSE]
    }
  }

  # split PIPs for molecular traits (e.g. introns) mapped to multiple genes
  if (any(duplicated(finemap_gene_res$id))) {
    loginfo("Split PIPs for molecular traits mapped to multiple genes")
    finemap_gene_res <- finemap_gene_res %>%
      group_by(.data$id) %>%
      mutate(susie_pip = ifelse(n() > 1, .data$susie_pip / n(), .data$susie_pip)) %>%
      ungroup()
  }

  # SNP results
  finemap_snp_res <- finemap_res[finemap_res$group=="SNP",]
  finemap_snp_res$gene_id <- NA
  finemap_snp_res$gene_name <- NA
  finemap_snp_res$gene_type <- "SNP"

  new_colnames <- unique(c("gene_id", "gene_name", "gene_type", colnames(finemap_res)))

  return(rbind(finemap_gene_res[,new_colnames],
               finemap_snp_res[,new_colnames]))
}


#' @title Add SNP and gene positions to cTWAS finemapping result.
#'
#' @param annotated_finemap_res a data frame of annotated cTWAS finemapping result
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param gene_annot a data frame of gene annotations, with required columns:
#' "gene_name", "chrom", "pos" (or "start" and "end").
#'
#' @param use_gene_pos if "pos" is not available in gene_annot,
#' use mid (midpoint), start or end positions to
#' represent gene positions.
#'
#' @param map_by column name to be mapped by (default: "gene_name").
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
add_pos_to_finemap_res <- function(annotated_finemap_res,
                                   snp_map,
                                   gene_annot,
                                   use_gene_pos = c("mid", "start", "end"),
                                   map_by = "gene_name"){

  loginfo("Annotating fine-mapping result ...")

  use_gene_pos <- match.arg(use_gene_pos)

  # gene results
  finemap_gene_res <- annotated_finemap_res[annotated_finemap_res$group!="SNP",]

  # Check required columns for gene_annot
  annot_cols <- c("gene_name", "start", "end")
  if (!all(annot_cols %in% colnames(gene_annot))){
    stop("gene_annot needs to contain the following columns: ",
         paste(annot_cols, collapse = " "))
  }
  gene_annot$chrom <- parse_number(as.character(gene_annot$chrom))
  gene_annot$start <- as.numeric(gene_annot$start)
  gene_annot$end <- as.numeric(gene_annot$end)

  # add gene positions
  loginfo("Add gene positions")
  # finemap_gene_res <- finemap_gene_res %>%
  #   left_join(gene_annot, by = map_by)
  gene_idx <- match(finemap_gene_res[,map_by], gene_annot[,map_by])
  finemap_gene_res$chrom <- gene_annot$chrom[gene_idx]

  # if "pos" is available in gene_annot, use "pos"
  if (!is.null(gene_annot$pos)){
    finemap_gene_res$pos <- gene_annot$pos[gene_idx]
  } else {
    if (use_gene_pos == "mid"){
      finemap_gene_res$pos <- round((gene_annot$start[gene_idx]+gene_annot$end[gene_idx])/2)
    } else if (use_gene_pos == "start") {
      finemap_gene_res$pos <- gene_annot$start[gene_idx]
    } else if (use_gene_pos == "end") {
      finemap_gene_res$pos <- gene_annot$end[gene_idx]
    }
  }

  # add SNP positions
  loginfo("Add SNP positions")
  finemap_snp_res <- annotated_finemap_res[annotated_finemap_res$group=="SNP",]
  snp_info <- as.data.frame(rbindlist(snp_map, idcol = "region_id"))
  snp_idx <- match(finemap_snp_res$id, snp_info$id)
  finemap_snp_res$chrom <- snp_info$chrom[snp_idx]
  finemap_snp_res$chrom <- parse_number(as.character(finemap_snp_res$chrom))
  finemap_snp_res$pos <- as.numeric(snp_info$pos[snp_idx])

  new_colnames <- unique(c("chrom", "pos", colnames(annotated_finemap_res)))

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
  gene_annot <- gene_annot[, c("chrom", "start", "end", "gene_id", "gene_name", "gene_type")]

  return(gene_annot)
}

#' @title Combines gene PIPs by context, type or group.
#'
#' @param finemap_res a data frame of annotated cTWAS finemapping result
#'
#' @param by sum gene PIPs by "context", "type", or "group".
#'
#' @param filter_protein_coding_genes If TRUE, keep protein coding genes only.
#'
#' @param filter_cs If TRUE, limits results in credible sets.
#'
#' @param missing_value set missing value (default: NA)
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
combine_gene_pips <- function(annotated_finemap_res,
                              by = c("context", "type", "group"),
                              gene_col = "gene_name",
                              missing_value = NA,
                              digits = 3){

  by <- match.arg(by)

  # Check to see if gene_name and gene_type are already in annotated_finemap_res
  if (!(gene_col %in% colnames(annotated_finemap_res))){
    stop("annotated_finemap_res needs to contain the column: ", gene_col,
         "\nPlease first run anno_finemap_res() to annotate finemap_res!")
  }

  # work with gene results below
  finemap_gene_res <- annotated_finemap_res[annotated_finemap_res$type!="SNP",]

  # combine PIPs
  combined_gene_pips <- aggregate(finemap_gene_res[,"susie_pip"],
                                  by = list(finemap_gene_res[,gene_col]),
                                  FUN = get_combined_pip)
  colnames(combined_gene_pips) <- c(gene_col, "combined_pip")

  if (by == "context"){
    # combine PIPs for each context
    contexts <- unique(finemap_gene_res$context)
    for (context in contexts){
      tmp_res <- finemap_gene_res[finemap_gene_res$context==context, c(gene_col, "susie_pip")]
      tmp_res <- aggregate(tmp_res[,"susie_pip"], by=list(tmp_res[,gene_col]), FUN=get_combined_pip)
      colnames(tmp_res) <- c(gene_col, paste0(context, "_pip"))
      combined_gene_pips <- combined_gene_pips %>% left_join(tmp_res, by = gene_col)
    }
  } else if (by == "type") {
    # combine PIPs for each type
    types <- unique(finemap_gene_res$type)
    for (type in types){
      tmp_res <- finemap_gene_res[finemap_gene_res$type==type, c(gene_col, "susie_pip")]
      tmp_res <- aggregate(tmp_res[,"susie_pip"], by=list(tmp_res[,gene_col]), FUN=get_combined_pip)
      colnames(tmp_res) <- c(gene_col, paste0(type, "_pip"))
      combined_gene_pips <- combined_gene_pips %>% left_join(tmp_res, by = gene_col)
    }
  } else if (by == "group") {
    # combine PIPs for each group
    groups <- unique(finemap_gene_res$group)
    for (group in groups){
      tmp_res <- finemap_gene_res[finemap_gene_res$group==group, c(gene_col, "susie_pip")]
      tmp_res <- aggregate(tmp_res[,"susie_pip"], by=list(tmp_res[,gene_col]), FUN=get_combined_pip)
      colnames(tmp_res) <- c(gene_col, paste0(group, "_pip"))
      combined_gene_pips <- combined_gene_pips %>% left_join(tmp_res, by = gene_col)
    }
  }

  if (!is.na(missing_value)) {
    combined_gene_pips[is.na(combined_gene_pips)] <- missing_value
  }

  # order by combined PIP
  combined_gene_pips <- combined_gene_pips[order(-combined_gene_pips$combined_pip),]

  # round gene PIPs
  combined_gene_pips[,-1] <- round(combined_gene_pips[, -1], digits)

  new_colnames <- c(setdiff(colnames(combined_gene_pips), "combined_pip"), "combined_pip")
  combined_gene_pips <- combined_gene_pips[, new_colnames]
  rownames(combined_gene_pips) <- NULL

  return(combined_gene_pips)
}

# combined gene PIP: 1 - \prod_k (1 - PIP_k).
# PIP_k is the PIP of the k-th molecular trait of a gene.
get_combined_pip <- function(pips){
  return(1 - prod(1 - pips))
}


#' Filter fine-mapping result
#'
#' @param annotated_finemap_res Annotated cTWAS fine-mapping result
#'
#' @param filter_protein_coding_genes If TRUE, keep protein coding genes only.
#'
#' @param filter_cs_genes If TRUE, limits gene results in credible sets.
#'
#' @param filter_cs_snps If TRUE, limits SNP results in credible sets.
#'
#' @return a data frame of filtered cTWAS fine-mapping result
#'
#' @export
filter_finemap_res <- function(annotated_finemap_res,
                               filter_protein_coding_genes = TRUE,
                               filter_cs_genes = TRUE,
                               filter_cs_snps = FALSE){

  finemap_gene_res <- annotated_finemap_res[annotated_finemap_res$type!="SNP",]

  # limit to protein coding genes
  if (filter_protein_coding_genes) {
    # Check to see if gene_name and gene_type are already in annotated_finemap_res
    if ("protein_coding" %in% colnames(annotated_finemap_res)){
      loginfo("Limit to protein coding genes")
      finemap_gene_res <- finemap_gene_res[finemap_gene_res$gene_type=="protein_coding",]
    } else {
      loginfo("No 'protein_coding' in 'gene_type'. Skipped filtering protein coding genes.")
    }
  }

  # limit genes in credible sets
  if (filter_cs_genes) {
    finemap_gene_res <- finemap_gene_res[finemap_gene_res$cs_index!=0,]
  }

  finemap_snp_res <- annotated_finemap_res[annotated_finemap_res$group=="SNP",]

  # limit SNPs in credible sets
  if (filter_cs_snps) {
    finemap_snp_res <- finemap_snp_res[finemap_snp_res$cs_index!=0,]
  }

  return(rbind(finemap_gene_res, finemap_snp_res))
}
