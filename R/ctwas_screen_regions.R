#' Screen regions with strong non-SNP PIPs to be used in finemapping step
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param weights a list of weights
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param L the number of effects for susie
#'
#' @param minvar minimum number of variables (snps and genes) in a region
#'
#' @param mingene minimum number of genes in a region
#'
#' @param min_nonSNP_PIP Regions with non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of screened_region_data, and a data frame of region ids and non-SNP PIPs
#'
#' @export
#'
screen_regions <- function(
    region_data,
    region_info,
    weights,
    group_prior = NULL,
    group_prior_var = NULL,
    L = 5,
    minvar = 2,
    mingene = 1,
    min_nonSNP_PIP = 0.5,
    ncore = 1,
    logfile = NULL,
    verbose = FALSE,
    ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Screening regions ...")

  # extract thin value from region_data
  thin <- unique(sapply(region_data, "[[", "thin"))
  if (length(thin) > 1) {
    thin <- min(thin)
    loginfo("thin has more than one value in region_data, use the minimum thin value.")
  }
  loginfo("thin = %s", thin)

  if (thin <= 0 || thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  # adjust to account for thin argument
  if (!is.null(group_prior)){
    group_prior["SNP"] <- group_prior["SNP"]/thin
  }

  # remove regions with fewer than minvar variables
  if (minvar > 0) {
    n.var <- sapply(region_data, function(x){length(x[["z"]])})
    drop.idx <- which(n.var < minvar)
    if (length(drop.idx) > 0){
      loginfo("Remove %d regions with number of variables < %d.", length(drop.idx), minvar)
      region_data[drop.idx] <- NULL
    }
  }

  # remove regions with fewer than mingene genes
  if (mingene > 0) {
    n.gid <- sapply(region_data, function(x){length(x[["gid"]])})
    drop.idx <- which(n.gid < mingene)
    if (length(drop.idx) > 0){
      loginfo("Remove %d regions with number of genes < %d.", length(drop.idx), mingene)
      region_data[drop.idx] <- NULL
    }
  }

  # run finemapping for all regions containing thinned SNPs
  loginfo("Run initial screening for %d regions ...", length(region_data))
  finemap_screening_res <- finemap_regions(region_data,
                                           region_info,
                                           weights,
                                           L = L,
                                           group_prior = group_prior,
                                           group_prior_var = group_prior_var,
                                           annotate_susie_result = FALSE,
                                           ncore = ncore,
                                           verbose = verbose,
                                           ...)

  # select regions based on total non-SNP PIP of the region
  loginfo("Computing non-SNP PIPs ...")
  region_nonSNP_PIP_df <- compute_region_nonSNP_PIPs(finemap_screening_res)

  loginfo("Select regions with non-SNP PIP >= %s", min_nonSNP_PIP)
  screened_region_ids <- region_nonSNP_PIP_df$region_id[region_nonSNP_PIP_df$nonSNP_PIP >= min_nonSNP_PIP]
  screened_region_data <- region_data[screened_region_ids]
  loginfo("Number of regions selected: %d", length(screened_region_data))

  return(list(screened_region_data = screened_region_data,
              region_nonSNP_PIP_df = region_nonSNP_PIP_df))
}

