#' @title Screens regions with strong non-SNP PIPs to be used in finemapping step
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_info a list of SNP info data frames for LD reference. Required when \code{use_LD = TRUE}.
#'
#' @param weights a list of weights for each gene. Required when \code{use_LD = TRUE}.
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
#' @param filter_L If TRUE, screening regions with L > 0
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
                           LD_info = NULL,
                           snp_info = NULL,
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

  estimated_L <- NULL
  # with-LD version: run finemapping for all regions containing thinned SNPs and estimate L for each region
  # no-LD version: set L = 1
  if (use_LD) {
    if (filter_L) {
      loginfo("Estimating L with uniform prior and screening regions with L > 0...")
      finemap_unif_prior_res <- finemap_regions(region_data,
                                                use_LD = TRUE,
                                                LD_info = LD_info,
                                                snp_info = snp_info,
                                                weights = weights,
                                                L = 5,
                                                include_cs_index = TRUE,
                                                LD_format = LD_format,
                                                LD_loader_fun = LD_loader_fun,
                                                ncore = ncore,
                                                verbose = verbose,
                                                ...)
      estimated_L <- get_L(finemap_unif_prior_res)
      screened_region_ids <- names(estimated_L[estimated_L > 0])
      screened_region_data <- region_data[screened_region_ids]
      loginfo("Number of regions selected with L > 0: %d", length(screened_region_data))
      L <- estimated_L[screened_region_ids]
    } else {
      loginfo("L = %d", L)
      screened_region_data <- region_data
    }

  } else {
    loginfo("Set L = 1 in no-LD version")
    L <- 1
    # use non-SNP PIP filtering in no-LD version
    filter_nonSNP_PIP <- TRUE
    screened_region_data <- region_data
  }

  nonSNP_PIPs <- NULL
  if (filter_nonSNP_PIP) {
    # if no-LD version: further filter regions by non-SNP PIP
    # filter regions with non-SNP PIPs for the regions with L > 0
    loginfo("Computing non-SNP PIPs and screening regions with non-SNP PIP >= %s", min_nonSNP_PIP)
    finemap_screening_res <- finemap_regions(screened_region_data,
                                             use_LD = use_LD,
                                             LD_info = LD_info,
                                             snp_info = snp_info,
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
    nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_screening_res)
    screened_region_ids <- names(nonSNP_PIPs[nonSNP_PIPs >= min_nonSNP_PIP])
    screened_region_data <- region_data[screened_region_ids]
    loginfo("Number of regions selected with non-SNP PIP >= %s: %d", min_nonSNP_PIP, length(screened_region_data))
    if (length(L) > 1) {
      L <- L[screened_region_ids]
    }
  }

  return(list(screened_region_data = screened_region_data,
              L = L,
              estimated_L = estimated_L,
              nonSNP_PIPs = nonSNP_PIPs))
}


#' @title Computes non-SNP PIPs for each region from finemapping result
#'
#' @param finemap_res a data frame of finemapping result
#'
#' @export
compute_region_nonSNP_PIPs <- function(finemap_res){
  region_ids <- unique(finemap_res$region_id)
  if (length(region_ids) == 0) {
    stop("no region_ids in finemap_res!")
  }
  nonSNP_PIPs <- sapply(region_ids, function(x){
    finemap_region_res <- finemap_res[finemap_res$region_id == x,]
    nonSNP_PIP <- sum(finemap_region_res$susie_pip[finemap_region_res$type != "SNP"])
    nonSNP_PIP[is.na(nonSNP_PIP)] <- 0 # 0 if nonSNP_PIP is NA
    nonSNP_PIP
  })
  names(nonSNP_PIPs) <- region_ids
  return(nonSNP_PIPs)
}

#' @title get L for each region from finemapping result
#'
#' @param finemap_res a data frame of finemapping result
#'
#' @export
get_L <- function(finemap_res){
  region_ids <- unique(finemap_res$region_id)
  if (length(region_ids) == 0) {
    stop("no region_ids in finemap_res!")
  }
  # get L for each region
  region_L <- sapply(region_ids, function(x){
    region_cs_index <- unique(finemap_res[finemap_res$region_id == x, "cs_index"])
    length(which(region_cs_index > 0))
  })
  names(region_L) <- region_ids
  return(region_L)
}
