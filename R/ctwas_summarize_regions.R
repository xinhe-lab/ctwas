
#' Estimate L for all regions by running finemapping with uniform prior
#'
#' @param region_data a list object indexing regions, variants and genes.
#'
#' @param LD_map a data frame with filenames of LD matrices for each of the regions.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param init_L upper bound of the number of causal signals
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a credible set.
#'
#' @param snps_only If TRUE, use only SNPs in the region data.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param verbose If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @importFrom logging loginfo
#'
#' @return a vector of estimated L for all regions
#'
#' @export
#'
estimate_region_L <- function(region_data,
                              LD_map,
                              weights,
                              init_L = 5,
                              min_abs_corr = 0.1,
                              snps_only = FALSE,
                              LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                              LD_loader_fun = NULL,
                              snpinfo_loader_fun = NULL,
                              ncore = 1,
                              verbose = FALSE,
                              ...) {

  # check inputs
  LD_format <- match.arg(LD_format)

  if (snps_only)
    loginfo("Use SNPs only when estimating L")

  # fine-mapping using uniform prior
  res <- finemap_regions(region_data,
                         LD_map = LD_map,
                         weights = weights,
                         L = init_L,
                         min_abs_corr = min_abs_corr,
                         include_cs = TRUE,
                         get_susie_alpha = FALSE,
                         snps_only = snps_only,
                         LD_format = LD_format,
                         LD_loader_fun = LD_loader_fun,
                         snpinfo_loader_fun = snpinfo_loader_fun,
                         ncore = ncore,
                         verbose = verbose,
                         ...)
  finemap_res <- res$finemap_res
  rm(res)
  all_estimated_L <- get_L(finemap_res)

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
    region_cs <- unique(na.omit(finemap_res[finemap_res$region_id == x, "cs"]))
    if (length(region_cs) > 0){
      n_cs <- length(unique(unlist(strsplit(region_cs, ","))))
    } else {
      n_cs <- 0
    }
    n_cs
  })
  names(region_L) <- region_ids
  return(region_L)
}


#' @title Computes non-SNP PIPs for all regions from finemapping result
#'
#' @param finemap_res a data frame of finemapping result
#'
#' @param filter_cs If TRUE, limits to credible sets.
#'
#' @return a vector of non-SNP PIPs for all regions
#'
#' @export
compute_region_nonSNP_PIPs <- function(finemap_res, filter_cs = TRUE){

  region_ids <- unique(finemap_res$region_id)
  if (length(region_ids) == 0) {
    stop("No region_ids in finemap_res to compute non-SNP PIPs!")
  }
  nonSNP_PIPs <- sapply(region_ids, function(x){
    finemap_region_res <- finemap_res[finemap_res$region_id == x,]
    if (filter_cs) {
      finemap_region_res <- finemap_region_res[!is.na(finemap_region_res$cs),,drop=FALSE]
    }
    nonSNP_PIP <- sum(finemap_region_res$susie_pip[finemap_region_res$group != "SNP"])
    nonSNP_PIP[is.na(nonSNP_PIP)] <- 0 # 0 if nonSNP_PIP is NA
    nonSNP_PIP
  })
  names(nonSNP_PIPs) <- region_ids
  return(nonSNP_PIPs)
}


# compute p(single effect)
compute_region_p_single_effect <- function(region_data, group_prior){

  if (!inherits(region_data,"list"))
    stop("'region_data' should be a list.")

  region_ids <- names(region_data)
  if (length(region_ids) == 0)
    stop("No region_ids in region_data!")

  groups <- unique(unlist(lapply(region_data, "[[", "groups")))
  groups_without_prior <- setdiff(groups, names(group_prior))
  if (length(groups_without_prior) > 0) {
    stop(paste("Missing group_prior for group:", groups_without_prior, "!"))
  }
  group_prior <- group_prior[names(group_prior) %in% groups]

  p_single_effect <- sapply(region_ids, function(region_id){
    gs_group <- extract_region_data(region_data, region_id)$gs_group
    group_size <- table(gs_group)[names(group_prior)]
    group_size[is.na(group_size)] <- 0
    p1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))
    p1
  })
  names(p_single_effect) <- region_ids

  return(p_single_effect)
}


