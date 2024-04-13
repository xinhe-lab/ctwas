#' Screen regions
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
#' @param mingene minimum number of genes in a region
#'
#' @param min_nonSNP_PIP Regions with non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a data frame of screened region ids and non-SNP PIPs
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

  # remove regions with fewer than mingene genes
  if (mingene > 0) {
    n.gid <- sapply(region_data, function(x){length(x[["gid"]])})
    drop.idx <- which(n.gid < mingene)
    loginfo("Remove %d regions with number of genes < %d.", length(drop.idx), mingene)
    region_data[drop.idx] <- NULL
  }

  # run finemapping for all regions containing thinned SNPs
  loginfo("Run initial screening for %d regions ...", length(region_data))
  finemap_res <- finemap_regions(region_data,
                                 region_info,
                                 weights,
                                 L = L,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 max_iter = max_iter,
                                 annotate_susie_result = FALSE,
                                 ncore = ncore,
                                 verbose = verbose)

  # select regions based on total non-SNP PIP of the region
  loginfo("Compute non-SNP PIPs")
  region_nonSNP_PIP_df <- compute_region_nonSNP_PIPs(finemap_res)
  loginfo("Select regions with non-SNP PIP >= %s", min_nonSNP_PIP)
  region_nonSNP_PIP_df <- region_nonSNP_PIP_df[region_nonSNP_PIP_df$nonSNP_PIP >= min_nonSNP_PIP, , drop=FALSE]
  loginfo("Number of regions selected: %d", nrow(region_nonSNP_PIP_df))

  return(region_nonSNP_PIP_df)
}

