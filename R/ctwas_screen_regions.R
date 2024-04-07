#' Screen regions
#'
#' @param regionlist a list object indexing regions, variants and genes.
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param weights a list of weights
#'
#' @param thin The proportion of SNPs to be used for parameter estimation and initial screening regions.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#'
#' @param mingene minimum number of genes in a region
#'
#' @param min_nonSNP_PIP Regions with non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param L the number of effects for susie
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a vector of screened region tags
#'
#' @export
#'
screen_regions <- function(
    regionlist,
    region_info,
    weights,
    thin = 1,
    group_prior = NULL,
    group_prior_var = NULL,
    L = 5,
    mingene = 1,
    min_nonSNP_PIP = 0.5,
    max_iter = 100,
    ncore = 1,
    logfile = NULL,
    verbose = FALSE,
    ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Screen regions ...")

  if (thin <= 0 || thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  if (thin != 1) {
    # adjust group_prior parameter to account for thin argument
    group_prior["SNP"] <- group_prior["SNP"] / thin
  }

  # remove regions with fewer than mingene genes
  if (mingene > 0) {
    n.gid <- sapply(regionlist, function(x){length(x[["gid"]])})
    drop.idx <- which(n.gid < mingene)
    loginfo("Remove %d regions with number of genes < %d.", length(drop.idx), mingene)
    regionlist[drop.idx] <- NULL
  }

  # run finemapping for all regions containing thinned SNPs
  loginfo("Run initial screening for %d regions ...", length(regionlist))
  finemap_res <- finemap_regions(regionlist,
                                 region_info,
                                 weights,
                                 L = L,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 annotate_susie_result = FALSE,
                                 ncore = ncore,
                                 verbose = verbose)

  # select regions based on total non-SNP PIP of the region
  loginfo("Select regions with non-SNP PIP >= %s", min_nonSNP_PIP)
  screened_region_tags <- select_highPIP_regions(finemap_res, min_nonSNP_PIP)
  loginfo("Number of regions selected: %d", length(screened_region_tags))

  return(screened_region_tags)
}

