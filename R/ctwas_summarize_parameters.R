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

  # estimated group prior inclusion probabilieis
  group_prior <- param$group_prior

  # estimated group prior variance
  group_prior_var <- param$group_prior_var

  # set group size
  group_size <- param$group_size
  group_size <- group_size[names(group_prior)]

  # estimated enrichment
  enrichment <- group_prior[names(group_prior) != "SNP"]/group_prior[names(group_prior) == "SNP"]

  # estimated group PVE
  group_pve <- group_prior_var*group_prior*group_size/gwas_n

  # total PVE
  total_pve <- sum(group_pve, na.rm = TRUE)

  # proportion of heritability
  prop_heritability <- group_pve/total_pve

  res <- list(group_size = group_size,
              group_prior = group_prior,
              group_prior_var = group_prior_var,
              enrichment = enrichment,
              group_pve = group_pve,
              total_pve = total_pve,
              prop_heritability = prop_heritability)

  return(res)
}
