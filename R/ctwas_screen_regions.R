#' Screen regions
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weight_list
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
    z_snp,
    z_gene,
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
    loginfo("Remove %d regions with fewer than %d genes.", length(drop.idx), mingene)
    regionlist[drop.idx] <- NULL
  }

  # run finemapping for all regions containing thinned SNPs
  loginfo("Run initial screening for %d regions ...", length(regionlist))
  finemap_res <- finemap_regions(z_snp,
                                 z_gene,
                                 regionlist,
                                 region_info,
                                 weights,
                                 L = L,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 ncore = ncore,
                                 verbose = verbose)

  # select regions based on total non-SNP PIP of the region
  loginfo("Select regions with non-SNP PIP >= %s", min_nonSNP_PIP)
  screened_region_tags <- NULL
  for (region_tag in names(regionlist)){
    finemap_region_res <- finemap_res[finemap_res$region_tag == region_tag,]
    nonSNP_PIP <- sum(finemap_region_res$susie_pip[finemap_region_res$type != "SNP"])
    nonSNP_PIP[is.na(nonSNP_PIP)] <- 0 # 0 if nonSNP_PIP is NA
    if (nonSNP_PIP >= min_nonSNP_PIP) {
      screened_region_tags <- c(screened_region_tags, region_tag)
    }
  }

  loginfo("Number of region tags that contain strong gene signals: %d", length(screened_region_tags))

  return(screened_region_tags)
}


