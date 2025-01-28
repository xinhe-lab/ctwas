#' Compute standard error and p-value for enrichment
#'
#' @param region_data a list of assembled region data.
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#'
#' @param group_prior_null a vector of prior inclusion probabilities for different groups under the null.
#'
#' @param group_prior_var a vector of prior variances for different groups.
#'
#' @return a list of estimated enrichment, S.E. and p-value.
#'
#' @export
#'
compute_enrichment_se <- function(region_data,
                                  group_prior,
                                  group_prior_null,
                                  group_prior_var){

  # estimated enrichment
  enrichment <- group_prior[names(group_prior) != "SNP"]/group_prior[names(group_prior) == "SNP"]

  loglik <- compute_loglik_ser(group_prior = group_prior,
                               group_prior_var = group_prior_var,
                               region_data = region_data,
                               null_method = "ctwas")

  # null enrichment (enrichment = 1)
  enrichment_null = 1
  loglik_null <- compute_loglik_ser(group_prior = group_prior_null,
                                    group_prior_var = group_prior_var,
                                    region_data = region_data,
                                    null_method = "ctwas")

  # compute S.E. using asymptotic test based on LRT
  # se = sqrt( theta^2 / (2 * (loglik - loglik_null)) )
  se_enrichment <- sqrt( enrichment^2 / (2 * (loglik - loglik_null)) )

  # p-value of enrichment
  z_enrichment <- (enrichment - enrichment_null) / se_enrichment
  p_enrichment <- 2*pnorm(abs(z_enrichment), lower.tail=FALSE)

  # # 95% confidence interval
  # ci_enrichment <- data.frame(lower.CI = enrichment - 1.96*se_enrichment,
  #                             upper.CI = enrichment + 1.96*se_enrichment)

  return(list("enrichment" = enrichment,
              "SE" = se_enrichment,
              "pval" = p_enrichment))
}




