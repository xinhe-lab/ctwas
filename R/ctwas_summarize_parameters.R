#' @title Summarizes estimated parameters
#'
#' @param param a list of parameter estimation result from \code{est_param}
#'
#' @param gwas_n the sample size of the GWAS summary statistics
#'
#' @return a list of summarized parameters
#'
#' @export
#'
summarize_param <- function(param, gwas_n){

  group_prior <- param$group_prior
  # estimated group prior (all iterations)
  group_prior_iters <- param$group_prior_iters

  group_prior_var <- param$group_prior_var
  # estimated group prior variance (all iterations)
  group_prior_var_iters <- param$group_prior_var_iters

  # set group size
  group_size <- param$group_size
  group_size <- group_size[rownames(group_prior_iters)]
  group_size <- as.numeric(group_size)
  names(group_size) <- rownames(group_prior_iters)

  # estimated group PVE (all iterations)
  group_pve_iters <- group_prior_var_iters*group_prior_iters*group_size/gwas_n
  group_pve <- group_pve_iters[,ncol(group_pve_iters)]

  # estimated enrichment of genes (all iterations)
  enrichment_iters <- t(sapply(rownames(group_prior_iters)[rownames(group_prior_iters)!="SNP"], function(x){
    group_prior_iters[rownames(group_prior_iters)==x,]/group_prior_iters[rownames(group_prior_iters)=="SNP"]}))
  enrichment <- enrichment_iters[,ncol(enrichment_iters)]

  res <- list(group_size = group_size,
              group_prior = group_prior,
              group_prior_var = group_prior_var,
              enrichment = enrichment,
              group_pve = group_pve,
              total_pve = sum(group_pve),
              prop_heritability = group_pve/sum(group_pve))

  return(res)
}
