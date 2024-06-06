#' Annotate cTWAS finemapping result with ensembl gene annotation database
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param snp_info a list or data frame of SNP info for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref".
#'
#' @param ens_db ensembl gene annotation database
#'
#' @param use_gene_pos use mid (midpoint), start or end positions to represent gene positions
#'
#' @importFrom data.table rbindlist
#' @importFrom ensembldb genes
#' @importFrom AnnotationFilter GeneIdFilter
#'
#' @return a data frame of cTWAS finemapping result including gene
#' names and updated positions
#'
#' @export
#' 
anno_finemap_res <- function(finemap_res, snp_info, ens_db,
                             use_gene_pos = c("mid", "start", "end")){

  use_gene_pos <- match.arg(use_gene_pos)

  # Check LD reference SNP info
  if (class(snp_info) == "list") {
    snp_info <- as.data.frame(data.table::rbindlist(snp_info, idcol = "region_id"))
  }
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("SNP info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  finemap_res <- cbind(chrom = NA, pos = NA, finemap_res, gene_id = NA, gene_name = NA)
  finemap_res <- finemap_res[, unique(colnames(finemap_res))]

  # extract gene IDs from finemapping result
  finemap_gene_res <- finemap_res[finemap_res$type!="SNP",]
  gene_id <- sapply(strsplit(finemap_gene_res$id, split = "[|]"), "[[", 1)
  gene_id <- sapply(strsplit(gene_id, split = "[.]"), "[[", 1)
  finemap_gene_res$gene_id <- gene_id

  # get gene position info
  gene_annot_gr <- genes(ens_db, filter = GeneIdFilter(finemap_gene_res$gene_id))
  gene_annot_gr$chrom <- seqnames(gene_annot_gr)
  gene_annot_gr$start_pos <- start(gene_annot_gr)
  gene_annot_gr$end_pos <- end(gene_annot_gr)
  if (use_gene_pos == "mid"){
    # use the midpoint as gene position
    gene_annot_gr$gene_pos <- round((gene_annot_gr$start_pos+gene_annot_gr$end_pos)/2)
  } else if (use_gene_pos == "start") {
    gene_annot_gr$gene_pos <- gene_annot_gr$start_pos
  } else if (use_gene_pos == "end") {
    gene_annot_gr$gene_pos <- gene_annot_gr$end_pos
  }

  # get gene_name and gene position
  gene_annot_df <- data.frame(gene_id = gene_annot_gr$gene_id,
                              gene_name = gene_annot_gr$gene_name,
                              gene_chrom = gene_annot_gr$chrom,
                              gene_pos = gene_annot_gr$gene_pos)

  # add gene name to finemapping result and update gene positions
  gene_annot_idx <- match(finemap_gene_res$gene_id, gene_annot_df$gene_id)
  finemap_gene_res$chrom <- gene_annot_df$gene_chrom[gene_annot_idx]
  finemap_gene_res$pos <- gene_annot_df$gene_pos[gene_annot_idx]
  finemap_gene_res$gene_name <- gene_annot_df$gene_name[gene_annot_idx]

  # add SNP positions
  finemap_SNP_res <- finemap_res[finemap_res$type=="SNP",]
  snp_idx <- match(finemap_SNP_res$id, snp_info$id)
  finemap_SNP_res$chrom <- snp_info$chrom[snp_idx]
  finemap_SNP_res$pos <- snp_info$pos[snp_idx]

  finemap_res <- rbind(finemap_gene_res, finemap_SNP_res)

  return(finemap_res)
}

#' Combine PIPs across contexts (tissues)
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param contexts a character vector of contexts to sum
#'
#' @return a data frame of combined gene PIPs and PIPs for each context
#' @export
sum_pip_across_contexts <- function(finemap_res, contexts){

  finemap_gene_res <- finemap_res[finemap_res$type!="SNP",]

  # combine PIPs across contexts (tissues)
  gene_pip_df <- aggregate(finemap_gene_res$susie_pip, by=list(finemap_gene_res$gene_id), FUN=sum)
  colnames(gene_pip_df) <- c("gene_id", "combined_pip")

  if (!is.null(finemap_gene_res$gene_name)){
    idx <- match(gene_pip_df$gene_id, finemap_gene_res$gene_id)
    gene_pip_df$gene_name <- finemap_gene_res$gene_name[idx]
    gene_pip_df <- gene_pip_df[, c("gene_id", "gene_name", "combined_pip")]
  }

  # PIPs for each context
  if (missing(contexts)){
    contexts <- unique(finemap_gene_res$context)
  }

  gene_pips_df <- matrix(NA, nrow=nrow(gene_pip_df), ncol=length(contexts))
  colnames(gene_pips_df) <- paste0(contexts, "_pip")
  for (i in 1:nrow(gene_pips_df)){
    gene_id <- gene_pip_df$gene_id[i]
    finemap_gene_res_subset <- finemap_gene_res[which(finemap_gene_res$gene_id==gene_id),]
    gene_contexts <- finemap_gene_res_subset$context
    gene_pips_df[i, paste0(gene_contexts, "_pip")] <- finemap_gene_res_subset$susie_pip
  }
  gene_pips_df <- cbind(gene_pip_df, gene_pips_df)

  # sort by combined PIP
  gene_pips_df <- gene_pips_df[order(-gene_pips_df$combined_pip),]

  return(gene_pips_df)
}
