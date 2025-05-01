
#' @title Screens regions with strong signals
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param min_var minimum number of variables (SNPs and genes) in a region.
#'
#' @param min_gene minimum number of genes in a region.
#'
#' @param min_nonSNP_PIP If screening by non-SNP PIPs,
#' regions with total non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param min_pval Keep regions with minimum p-values from z_snp and z_gene < \code{min_pval}.
#'
#' @param null_method Method to compute null model, options: "ctwas", "susie" or "none".
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1). Only used when \code{null_method = "susie"}.
#'
#' @param ncore The number of cores used to parallelize susie over regions.
#'
#' @param logfile The log filename. If NULL, print log info on screen.
##'
#' @param verbose If TRUE, print detail messages.
#'
#' @return a list, containing a data frame of selected region data,
#' and a data frame of screening summary for all regions.
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#'
#' @export
#'
screen_regions <- function(region_data,
                           group_prior = NULL,
                           group_prior_var = NULL,
                           min_var = 2,
                           min_gene = 1,
                           min_nonSNP_PIP = 0.5,
                           min_pval = 5e-8,
                           null_method = c("ctwas", "susie", "none"),
                           null_weight = NULL,
                           ncore = 1,
                           logfile = NULL,
                           verbose = FALSE){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Screening regions ...")

  null_method <- match.arg(null_method)

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (anyDuplicated(names(region_data)))
    logwarn("Duplicated names of region_data found! Please use unique names for region_data!")

  # use all SNPs (thin = 1) for screening
  thin <- min(sapply(region_data, "[[", "thin"))
  if (thin != 1)
    stop("thin != 1, please run expand_region_data() first to include all SNPs!")

  # adjust group_prior to account for thin argument
  if (!is.null(group_prior)){
    group_prior["SNP"] <- group_prior["SNP"]/thin
  }

  # create a data frame for screening summary
  region_ids <- names(region_data)
  screen_summary <- summarize_region_signals(region_data)
  n_gids <- screen_summary$n_gids
  n_sids <- screen_summary$n_sids
  screen_summary$nonSNP_PIP <- NA

  # skip regions with fewer than min_var variables
  skipped_region_ids <- NULL
  if (min_var > 0) {
    min_var_region_ids <- region_ids[(n_sids + n_gids) < min_var]
    if (length(min_var_region_ids) > 0){
      loginfo("Skip %d regions with number of variables < %d.", length(min_var_region_ids), min_var)
      skipped_region_ids <- c(skipped_region_ids, min_var_region_ids)
    }
  }

  # skip regions with fewer than min_gene genes
  if (min_gene > 0) {
    min_gene_region_ids <- region_ids[n_gids < min_gene]
    if (length(min_gene_region_ids) > 0){
      loginfo("Skip %d regions with number of genes < %d.", length(min_gene_region_ids), min_gene)
      skipped_region_ids <- c(skipped_region_ids, min_gene_region_ids)
    }
  }

  sigP_region_ids <- screen_summary$region_id[which(screen_summary$min_gene_p < min_pval | screen_summary$min_snp_p < min_pval)]
  sigP_region_ids <- setdiff(sigP_region_ids, skipped_region_ids)
  loginfo("Selected %d regions with minimum p-value < %s.", length(sigP_region_ids), min_pval)

  # run finemapping for all regions without LD (L=1) using the SER model
  # select regions with total non-SNP PIPs > 0.5
  region_ids_to_screen <- setdiff(region_ids, c(skipped_region_ids, sigP_region_ids))
  loginfo("Screening %d regions ...", length(region_ids_to_screen))
  finemap_screening_res <- finemap_regions_ser(region_data[region_ids_to_screen],
                                               group_prior = group_prior,
                                               group_prior_var = group_prior_var,
                                               null_method = null_method,
                                               null_weight = null_weight,
                                               ncore = ncore,
                                               verbose = verbose)
  # select regions based on total non-SNP PIPs
  all_nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_screening_res, filter_cs = FALSE)
  idx <- match(names(all_nonSNP_PIPs), screen_summary$region_id)
  screen_summary$nonSNP_PIP[idx] <- all_nonSNP_PIPs
  screened_region_ids <- screen_summary$region_id[which(screen_summary$nonSNP_PIP > min_nonSNP_PIP)]
  loginfo("Selected %d regions with non-SNP PIP > %s", length(screened_region_ids), min_nonSNP_PIP)
  all_selected_region_ids <- sort(unique(c(screened_region_ids, sigP_region_ids)))
  screened_region_data <- region_data[all_selected_region_ids]
  loginfo("Selected %d regions in total", length(screened_region_data))

  return(list("screened_region_data" = screened_region_data,
              "screen_summary" = screen_summary))
}

