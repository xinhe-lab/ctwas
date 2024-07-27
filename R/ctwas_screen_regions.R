#' @title Screens regions with strong non-SNP PIPs to be used in finemapping step
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_map a data frame with filenames of LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of weights for each gene. Required when \code{use_LD = TRUE}.
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param L the number of effects for susie.
#'
#' @param minvar minimum number of variables (snps and genes) in a region
#'
#' @param mingene minimum number of genes in a region
#'
#' @param filter_L If TRUE, screening regions with estimated L > 0.
#'
#' @param filter_nonSNP_PIP If TRUE, screening regions with total non-SNP PIP >= \code{min_nonSNP_PIP}
#'
#' @param min_nonSNP_PIP If screening by non-SNP PIPs,
#' regions with total non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
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
#' @return a list of screened_region_data, and
#' non-SNP PIPs or estimated L for each region
#'
#' @export
#'
screen_regions <- function(region_data,
                           use_LD = TRUE,
                           LD_map = NULL,
                           snp_map = NULL,
                           weights = NULL,
                           group_prior = NULL,
                           group_prior_var = NULL,
                           L = 5,
                           minvar = 2,
                           mingene = 1,
                           filter_L = TRUE,
                           filter_nonSNP_PIP = FALSE,
                           min_nonSNP_PIP = 0.5,
                           LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                           LD_loader_fun,
                           ncore = 1,
                           logfile = NULL,
                           verbose = FALSE,
                           ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Screening regions ...")

  # check input
  LD_format <- match.arg(LD_format)

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

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (use_LD) {
    if (is.null(LD_map) || is.null(snp_map) || is.null(weights))
      stop("LD_map, snp_map and weights are required when use_LD = TRUE")

    if (!inherits(snp_map,"list"))
      stop("'snp_map' should be a list.")

    if (!inherits(LD_map,"data.frame"))
      stop("'LD_map' should be a data frame")

    if (!inherits(weights,"list"))
      stop("'weights' should be a list.")

    if (any(sapply(weights, is.null)))
      stop("weights contain NULL, remove empty weights!")
  }

  # adjust to account for thin argument
  if (!is.null(group_prior)){
    group_prior["SNP"] <- group_prior["SNP"]/thin
  }

  # remove regions with fewer than minvar variables
  if (minvar > 0) {
    n.var <- sapply(region_data, function(x){length(x$gid) + length(x$sid)})
    drop.idx <- which(n.var < minvar)
    if (length(drop.idx) > 0){
      loginfo("Remove %d regions with number of variables < %d.", length(drop.idx), minvar)
      region_data[drop.idx] <- NULL
    }
  }

  # remove regions with fewer than mingene genes
  if (mingene > 0) {
    n.gid <- sapply(region_data, function(x){length(x$gid)})
    drop.idx <- which(n.gid < mingene)
    if (length(drop.idx) > 0){
      loginfo("Remove %d regions with number of genes < %d.", length(drop.idx), mingene)
      region_data[drop.idx] <- NULL
    }
  }

  # no-LD version: set L = 1
  if (!use_LD) {
    L <- 1
    filter_nonSNP_PIP <- TRUE
    filter_L <- FALSE
    screened_region_data <- region_data
    all_estimated_L <- NULL
  } else {
    # with-LD version: run finemapping for all regions containing thinned SNPs and estimate L for each region
    if (filter_L) {
      loginfo("Estimating L with uniform prior ...")
      all_estimated_L <- estimate_region_L(region_data = region_data,
                                           LD_map = LD_map,
                                           snp_map = snp_map,
                                           weights = weights,
                                           LD_format = LD_format,
                                           LD_loader_fun = LD_loader_fun,
                                           ncore = ncore,
                                           verbose = verbose,
                                           ...)
      screened_region_ids <- names(all_estimated_L[all_estimated_L >= 1])
      screened_region_data <- region_data[screened_region_ids]
      loginfo("Selected %d regions with L >= 1", length(screened_region_data))
      L <- all_estimated_L[screened_region_ids]
    } else {
      L <- L
      loginfo("Set L = %d", L)
      screened_region_data <- region_data
      all_estimated_L <- NULL
    }
  }

  # filter regions with non-SNP PIPs for the regions
  if (filter_nonSNP_PIP) {
    loginfo("Computing non-SNP PIPs ...")
    finemap_screening_res <- finemap_regions(screened_region_data,
                                             use_LD = use_LD,
                                             LD_map = LD_map,
                                             snp_map = snp_map,
                                             weights = weights,
                                             L = L,
                                             group_prior = group_prior,
                                             group_prior_var = group_prior_var,
                                             include_cs_index = FALSE,
                                             LD_format = LD_format,
                                             LD_loader_fun = LD_loader_fun,
                                             ncore = ncore,
                                             verbose = verbose,
                                             ...)
    # select regions based on total non-SNP PIP of the region
    all_nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_screening_res)
    screened_region_ids <- names(all_nonSNP_PIPs[all_nonSNP_PIPs >= min_nonSNP_PIP])
    screened_region_data <- region_data[screened_region_ids]
    loginfo("Selected %d regions with non-SNP PIP >= %s", length(screened_region_data), min_nonSNP_PIP)
    if (length(L) > 1) {
      L <- L[screened_region_ids]
    }
  } else{
    all_nonSNP_PIPs <- NULL
  }

  return(list(screened_region_data = screened_region_data,
              screened_region_ids = screened_region_ids,
              L = L,
              all_estimated_L = all_estimated_L,
              all_nonSNP_PIPs = all_nonSNP_PIPs))
}



#' Estimate L for each region by running finemapping with uniform prior
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param LD_map a data frame with filenames of LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of weights for each gene. Required when \code{use_LD = TRUE}.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging loginfo
#'
#' @return estimated L for each region
#'
#' @export
#'
estimate_region_L <- function(region_data,
                              LD_map,
                              snp_map,
                              weights,
                              LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                              LD_loader_fun,
                              ncore = 1,
                              verbose = FALSE,
                              ...) {

  LD_format <- match.arg(LD_format)

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (!inherits(snp_map,"list"))
    stop("'snp_map' should be a list.")

  if (!inherits(LD_map,"data.frame"))
    stop("'LD_map' should be a data frame")

  if (!inherits(weights,"list"))
    stop("'weights' should be a list.")

  if (any(sapply(weights, is.null)))
    stop("weights contain NULL, remove empty weights!")

  finemap_unif_prior_res <- finemap_regions(region_data,
                                            use_LD = TRUE,
                                            LD_map = LD_map,
                                            snp_map = snp_map,
                                            weights = weights,
                                            L = 5,
                                            include_cs_index = TRUE,
                                            LD_format = LD_format,
                                            LD_loader_fun = LD_loader_fun,
                                            ncore = ncore,
                                            verbose = verbose,
                                            ...)
  all_estimated_L <- get_L(finemap_unif_prior_res)

  return(all_estimated_L)
}

# get L for each region from finemapping result
get_L <- function(finemap_res){
  region_ids <- unique(finemap_res$region_id)
  if (length(region_ids) == 0) {
    stop("No region_ids in finemap_res!")
  }
  # get L for each region
  region_L <- sapply(region_ids, function(x){
    region_cs_index <- unique(finemap_res[finemap_res$region_id == x, "cs_index"])
    length(which(region_cs_index > 0))
  })
  names(region_L) <- region_ids
  return(region_L)
}


#' @title Computes non-SNP PIPs for each region from finemapping result
#'
#' @param finemap_res a data frame of finemapping result
#'
#' @export
compute_region_nonSNP_PIPs <- function(finemap_res){
  region_ids <- unique(finemap_res$region_id)
  if (length(region_ids) == 0) {
    stop("No region_ids in finemap_res!")
  }
  nonSNP_PIPs <- sapply(region_ids, function(x){
    finemap_region_res <- finemap_res[finemap_res$region_id == x,]
    nonSNP_PIP <- sum(finemap_region_res$susie_pip[finemap_region_res$group != "SNP"])
    nonSNP_PIP[is.na(nonSNP_PIP)] <- 0 # 0 if nonSNP_PIP is NA
    nonSNP_PIP
  })
  names(nonSNP_PIPs) <- region_ids
  return(nonSNP_PIPs)
}
