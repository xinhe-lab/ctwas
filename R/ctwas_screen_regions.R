#' @title Screens regions with strong signals to be used in fine-mapping step.
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for each of the regions.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param L the number of effects for susie.
#'
#' @param min_var minimum number of variables (SNPs and genes) in a region.
#'
#' @param min_gene minimum number of genes in a region.
#'
#' @param filter_L If TRUE, screening regions with estimated credible sets
#' (L) = \code{min_L}.
#'
#' @param filter_nonSNP_PIP If TRUE, screening regions with
#' total non-SNP PIP >= \code{min_nonSNP_PIP}.
#'
#' @param min_L If screening by nubmer of credible sets (L),
#' regions with L >= \code{min_L}
#' will be selected to run finemapping using full SNPs.
#'
#' @param min_nonSNP_PIP If screening by non-SNP PIPs,
#' regions with total non-SNP PIP >= \code{min_nonSNP_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a credible set.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix
#' when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore The number of cores used to parallelize susie over regions.
#'
#' @param logfile The log filename. If NULL, will print log info on screen.
#'
#' @param verbose If TRUE, print detail messages.
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a list, containing a data frame of selected region data,
#' estimated L for selected regions,
#' and a data frame of screening summary for all regions.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @export
#'
screen_regions <- function(region_data,
                           LD_map,
                           weights,
                           group_prior = NULL,
                           group_prior_var = NULL,
                           L = 5,
                           min_var = 2,
                           min_gene = 1,
                           filter_L = TRUE,
                           filter_nonSNP_PIP = FALSE,
                           min_L = 1,
                           min_nonSNP_PIP = 0.5,
                           min_abs_corr = 0.1,
                           LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                           LD_loader_fun = NULL,
                           snpinfo_loader_fun = NULL,
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
    loginfo("thin has more than one value in region_data, use the minimum thin value")
  }
  loginfo("thin = %s", thin)

  if (thin <= 0 || thin > 1)
    stop("thin value needs to be in (0,1]")

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  if (!inherits(weights,"list"))
    stop("'weights' should be a list!")

  if (any(sapply(weights, is.null)))
    stop("'weights' contain NULL, remove empty weights!")

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

  # skip regions with fewer than min_var variables
  if (min_var > 0) {
    skip_region_ids <- region_ids[(n_sids + n_gids) < min_var]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of variables < %d", length(skip_region_ids), min_var)
      region_data[skip_region_ids] <- NULL
    }
  }

  # skip regions with fewer than min_gene genes
  if (min_gene > 0) {
    skip_region_ids <- region_ids[n_gids < min_gene]
    if (length(skip_region_ids) > 0){
      loginfo("Remove %d regions with number of genes < %d", length(skip_region_ids), min_gene)
      region_data[skip_region_ids] <- NULL
    }
  }

  # run finemapping for all regions with thinned SNPs using uniform prior
  # estimate L for each region
  # select regions with L > 0
  if (filter_L) {
    loginfo("Estimating L ...")
    all_estimated_L <- estimate_region_L(region_data = region_data,
                                         LD_map = LD_map,
                                         weights = weights,
                                         init_L = L,
                                         min_abs_corr = min_abs_corr,
                                         snps_only = FALSE,
                                         LD_format = LD_format,
                                         LD_loader_fun = LD_loader_fun,
                                         snpinfo_loader_fun = snpinfo_loader_fun,
                                         ncore = ncore,
                                         verbose = verbose,
                                         ...)
    screened_region_ids <- names(all_estimated_L[all_estimated_L >= min_L])
    screened_region_data <- region_data[screened_region_ids]
    loginfo("Selected %d regions with L >= %d", length(screened_region_data), min_L)
    screened_region_L <- all_estimated_L[screened_region_ids]
    idx <- match(names(all_estimated_L), screen_summary$region_id)
    screen_summary$L[idx] <- all_estimated_L
  } else {
    loginfo("L = %d", L)
    screened_region_data <- region_data
    screened_region_L <- L
    screen_summary$L <- L
  }

  # run finemapping for all regions with thinned SNPs using estimated priors
  # select regions with total non-SNP PIPs > 0.5
  if (filter_nonSNP_PIP) {
    loginfo("Computing non-SNP PIPs ...")
    res <- finemap_regions(screened_region_data,
                           LD_map = LD_map,
                           weights = weights,
                           L = L,
                           group_prior = group_prior,
                           group_prior_var = group_prior_var,
                           include_cs = FALSE,
                           get_susie_alpha = FALSE,
                           LD_format = LD_format,
                           LD_loader_fun = LD_loader_fun,
                           snpinfo_loader_fun = snpinfo_loader_fun,
                           ncore = ncore,
                           verbose = verbose,
                           ...)
    finemap_screening_res <- res$finemap_res
    rm(res)
    # select regions based on total non-SNP PIPs
    all_nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_screening_res, filter_cs = FALSE)
    screened_region_ids <- names(all_nonSNP_PIPs[all_nonSNP_PIPs >= min_nonSNP_PIP])
    screened_region_data <- region_data[screened_region_ids]
    loginfo("Selected %d regions with non-SNP PIP >= %s", length(screened_region_data), min_nonSNP_PIP)

    if (length(screened_region_L) > 1) {
      screened_region_L <- screened_region_L[screened_region_ids]
    }

    idx <- match(names(all_nonSNP_PIPs), screen_summary$region_id)
    screen_summary$nonSNP_PIP[idx] <- all_nonSNP_PIPs
  }

  # if min_L = 0, set regions with L = 0 to 1
  if (min_L == 0) {
    screened_region_L[screened_region_L == 0] <- 1
  }

  rownames(screen_summary) <- NULL

  return(list("screened_region_data" = screened_region_data,
              "screened_region_L" = screened_region_L,
              "screen_summary" = screen_summary))
}


#' @title Screens regions with strong signals without using LD.
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
#' @param ncore The number of cores used to parallelize susie over regions.
#'
#' @param logfile The log filename. If NULL, print log info on screen.
##'
#' @param verbose If TRUE, print detail messages.
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a list, containing a data frame of selected region data,
#' and a data frame of screening summary for all regions.
#'
#' @importFrom logging addHandler loginfo logwarn writeToFile
#'
#' @export
#'
screen_regions_noLD <- function(region_data,
                                group_prior = NULL,
                                group_prior_var = NULL,
                                min_var = 2,
                                min_gene = 1,
                                min_nonSNP_PIP = 0.5,
                                ncore = 1,
                                logfile = NULL,
                                verbose = FALSE,
                                ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Screening regions without LD ...")

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

  if (anyDuplicated(names(region_data)))
    logwarn("Duplicated names of region_data found! Please use unique names for region_data!")

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
                               L = 1,
                               nonSNP_PIP = NA)

  # skip regions with fewer than min_var variables
  if (min_var > 0) {
    skip_region_ids <- region_ids[(n_sids + n_gids) < min_var]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of variables < %d.", length(skip_region_ids), min_var)
      region_data <- region_data[!names(region_data) %in% skip_region_ids]
    }
  }

  # skip regions with fewer than min_gene genes
  if (min_gene > 0) {
    skip_region_ids <- region_ids[n_gids < min_gene]
    if (length(skip_region_ids) > 0){
      loginfo("Skip %d regions with number of genes < %d.", length(skip_region_ids), min_gene)
      region_data <- region_data[!names(region_data) %in% skip_region_ids]
    }
  }

  # run finemapping for all regions with thinned SNPs using estimated priors
  # select regions with total non-SNP PIPs > 0.5
  loginfo("Computing non-SNP PIPs ...")
  res <- finemap_regions_noLD(region_data,
                              group_prior = group_prior,
                              group_prior_var = group_prior_var,
                              get_susie_alpha = FALSE,
                              ncore = ncore,
                              verbose = verbose,
                              ...)
  finemap_screening_res <- res$finemap_res
  rm(res)
  # select regions based on total non-SNP PIPs
  all_nonSNP_PIPs <- compute_region_nonSNP_PIPs(finemap_screening_res, filter_cs = FALSE)
  screened_region_ids <- names(all_nonSNP_PIPs[all_nonSNP_PIPs >= min_nonSNP_PIP])
  screened_region_data <- region_data[screened_region_ids]
  loginfo("Selected %d regions with non-SNP PIP >= %s", length(screened_region_data), min_nonSNP_PIP)
  idx <- match(names(all_nonSNP_PIPs), screen_summary$region_id)
  screen_summary$nonSNP_PIP[idx] <- all_nonSNP_PIPs

  rownames(screen_summary) <- NULL

  return(list("screened_region_data" = screened_region_data,
              "screen_summary" = screen_summary))
}


