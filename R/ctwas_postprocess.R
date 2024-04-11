#' run cTWAS finemapping for multiple regions with L = 1
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param regionlist a list object indexing regions, variants and genes.
#'
#' @param region_id a character string of region id to be finemapped
#'
#' @param weight_list a list of weights for each gene
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @importFrom logging loginfo
#'
#' @return finemapping results.
#'
#' @export
#'
finemap_regions_L1 <- function(z_snp,
                               z_gene,
                               gene_info,
                               regionlist,
                               region_ids,
                               weight_list = NULL,
                               group_prior = NULL,
                               group_prior_var = NULL,
                               use_null_weight = TRUE,
                               max_iter = 1,
                               logfile = NULL,
                               ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  if (length(region_ids) == 0) {
    loginfo("No regions included!")
    finemap_res <- NULL
  }else{
    # select and assemble regionlist for rerunning finemapping
    loginfo('Finemmapping with L = 1 for %d regions...', length(region_ids))

    # Rerun finemapping with L = 1
    finemap_res <- list()
    for (region_id in region_ids) {
      finemap_res[[region_id]] <- finemap_region(z_snp = z_snp,
                                                  z_gene = z_gene,
                                                  gene_info = gene_info,
                                                  regionlist = regionlist,
                                                  region_id = region_id,
                                                  weight_list = weight_list,
                                                  L = 1,
                                                  group_prior = group_prior,
                                                  group_prior_var = group_prior_var,
                                                  use_null_weight = use_null_weight,
                                                  ...)
    }
    finemap_res <- do.call(rbind, finemap_res)
  }

  return(finemap_res)
}

#' Get regions with problematic high PIP SNPs or genes
#'
#' @param ctwas_res a data frame of cTWAS finemapping result
#' @param weight a string, weight filename used in cTWAS
#' @param problematic_snps a character vector of problematic SNP rsIDs
#' @param pip_thresh Minimum PIP value to select regions
#'
#' @return a character vector of region ids with problematic high PIP SNPs or genes
#'
#' @importFrom logging loginfo
#'
#' @export
get_problematic_regions <- function(ctwas_res, weight, problematic_snps, pip_thresh = 0.5){

  if (length(problematic_snps) == 0) {
    loginfo('No problematic SNPs')
    problematic_region_ids <- NULL
  }else{
    loginfo('Number of problematic SNPs: %d', length(problematic_snps))

    # read the PredictDB weights
    stopifnot(file.exists(weight))
    sqlite <- RSQLite::dbDriver("SQLite")
    db <- RSQLite::dbConnect(sqlite, weight)
    query <- function(...) RSQLite::dbGetQuery(db, ...)
    weight_table <- query("select * from weights")
    # load gene information from PredictDB weights
    gene_info <- query("select gene, genename, gene_type from extra")
    RSQLite::dbDisconnect(db)

    # find regions with high PIP results
    ctwas_highpip_res <- ctwas_res[ctwas_res$susie_pip > pip_thresh, ]

    # find regions with high PIP variants (PIP > 0.5) that are problematic
    ctwas_highpip_snp_res <- ctwas_highpip_res[ctwas_highpip_res$type == "SNP", ]
    problematic_highpip_snps <- intersect(ctwas_highpip_snp_res$id, problematic_snps)
    loginfo('Number of problematic high PIP SNPs: %d', length(problematic_highpip_snps))

    # find high PIP genes with problematic variants in its weights
    ctwas_highpip_gene_res <- ctwas_highpip_res[ctwas_highpip_res$type != "SNP", ]
    ctwas_highpip_gene_weight_table <- weight_table[weight_table$gene %in% ctwas_highpip_gene_res$id, ]
    problematic_highpip_genes <- ctwas_highpip_gene_weight_table$gene[which(ctwas_highpip_gene_weight_table$rsid %in% problematic_snps)]
    loginfo('Number of problematic high PIP genes: %d', length(problematic_highpip_genes))

    # get problematic high PIP regions
    problematic_ids <- c(problematic_highpip_snps, problematic_highpip_genes)
    if (length(problematic_ids) > 0) {
      problematic_region_ids <- unique(ctwas_res[ctwas_res$id %in% problematic_ids,"region_id"])
      loginfo('Number of problematic regions: %d', length(problematic_region_ids))
    }else{
      loginfo('No problematic regions found')
    }
  }

  # return problematic region ids
  return(problematic_region_ids)
}

