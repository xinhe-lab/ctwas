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
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
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
#' @param expand If TRUE, expand screened region data with all SNPs
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
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
                           LD_map,
                           snp_map,
                           weights,
                           z_snp,
                           group_prior = NULL,
                           group_prior_var = NULL,
                           L = 5,
                           minvar = 2,
                           mingene = 1,
                           filter_L = TRUE,
                           filter_nonSNP_PIP = FALSE,
                           min_nonSNP_PIP = 0.5,
                           expand = TRUE,
                           maxSNP = Inf,
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

  if (thin <= 0 || thin > 1)
    stop("thin value needs to be in (0,1]")

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (use_LD) {
    if (missing(LD_map) || missing(snp_map) || missing(weights))
      stop("'LD_map', 'snp_map' and 'weights' are required when use_LD = TRUE")

    if (!inherits(snp_map,"list"))
      stop("'snp_map' should be a list.")

    if (!inherits(LD_map,"data.frame"))
      stop("'LD_map' should be a data frame")

    if (!inherits(weights,"list"))
      stop("'weights' should be a list.")

    if (any(sapply(weights, is.null)))
      stop("'weights' contain NULL, remove empty weights!")
  }

  if (expand) {
    if (missing(snp_map) || missing(z_snp))
      stop("'snp_map' and 'z_snp' are required when expand = TRUE")
  }

  # adjust group_prior to account for thin argument
  if (!is.null(group_prior)){
    group_prior["SNP"] <- group_prior["SNP"]/thin
  }

  # create a data frame for screening summary
  region_ids <- names(region_data)
  n_gids <- sapply(region_data, function(x){length(x$gid)})
  n_sids <- sapply(region_data, function(x){length(x$sid)})
  screen_summary <- data.frame(region_id = region_ids,
                               n_gids = n_gids,
                               n_sids = n_sids,
                               L = NA,
                               nonSNP_PIP = NA)

  # skip regions with fewer than mingene genes
  if (mingene > 0) {
    skip_region_ids <- region_ids[n_gids < mingene]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of genes < %d.", length(skip_region_ids), mingene)
      region_data[skip_region_ids] <- NULL
    }
  }

  # skip regions with fewer than minvar variables
  if (minvar > 0) {
    skip_region_ids <- region_ids[(n_gids + n_sids) < minvar]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of variables < %d.", length(skip_region_ids), minvar)
      region_data[skip_region_ids] <- NULL
    }
  }

  # no-LD version: set L = 1
  if (!use_LD) {
    L <- 1
    screen_summary$L <- 1
    filter_nonSNP_PIP <- TRUE
    filter_L <- FALSE
    screened_region_data <- region_data
  } else {
    # with-LD version: run finemapping for all regions containing thinned SNPs and estimate L for each region
    if (filter_L) {
      loginfo("Estimating L ...")
      all_estimated_L <- estimate_region_L(region_data = region_data,
                                           LD_map = LD_map,
                                           snp_map = snp_map,
                                           weights = weights,
                                           init_L = L,
                                           LD_format = LD_format,
                                           LD_loader_fun = LD_loader_fun,
                                           ncore = ncore,
                                           verbose = verbose,
                                           ...)
      screened_region_ids <- names(all_estimated_L[all_estimated_L > 0])
      screened_region_data <- region_data[screened_region_ids]
      loginfo("Selected %d regions with L > 0", length(screened_region_data))
      L <- all_estimated_L[screened_region_ids]
      idx <- match(names(all_estimated_L), screen_summary$region_id)
      screen_summary$L[idx] <- all_estimated_L
    } else {
      loginfo("L = %d", L)
      screen_summary$L <- L
      screened_region_data <- region_data
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
    # select regions based on total non-SNP PIPs
    all_nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_screening_res)
    screened_region_ids <- names(all_nonSNP_PIPs[all_nonSNP_PIPs >= min_nonSNP_PIP])
    screened_region_data <- region_data[screened_region_ids]
    loginfo("Selected %d regions with non-SNP PIP >= %s", length(screened_region_data), min_nonSNP_PIP)
    idx <- match(names(all_nonSNP_PIPs), screen_summary$region_id)
    screen_summary$nonSNP_PIP[idx] <- all_nonSNP_PIPs

    if (length(L) > 1) {
      L <- L[screened_region_ids]
    }
  }

  if (expand) {
    if (thin < 1) {
      screened_region_data <- expand_region_data(screened_region_data,
                                                 snp_map = snp_map,
                                                 z_snp = z_snp,
                                                 maxSNP = maxSNP,
                                                 ncore = ncore)
    }
  }

  return(list("screened_region_data" = screened_region_data,
              "L" = L,
              "screen_summary" = screen_summary))
}



#' Estimate L for each region by running finemapping with uniform prior
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
                              init_L = 5,
                              LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                              LD_loader_fun,
                              ncore = 1,
                              verbose = FALSE,
                              ...) {

  # check inputs
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

  # fine-mapping using uniform prior
  finemap_unif_prior_res <- finemap_regions(region_data,
                                            use_LD = TRUE,
                                            LD_map = LD_map,
                                            snp_map = snp_map,
                                            weights = weights,
                                            L = init_L,
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
