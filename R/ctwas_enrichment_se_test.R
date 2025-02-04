#' Compute standard error and p-value for enrichment using
#' the likelihood ratio test from EM results
#'
#' @param region_data a list of assembled region data.
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#'
#' @param group_prior_var a vector of prior variances for different groups.
#'
#' @param niter the number of iterations of the E-M algorithm to perform
#'
#' @param null_method Method to compute null model, options: "ctwas", "susie" or "none".
#'
#' @param ncore The number of cores used to parallelize over regions
#'
#' @param log_enrichment If TRUE, return enrichment result in log scale.
#'
#' @return Estimated enrichment, S.E. and p-value.
#'
#' @export
#'
get_enrichment_se_LRT <- function(region_data,
                                  group_prior,
                                  group_prior_var,
                                  niter = 20,
                                  null_method = c("ctwas", "susie", "none"),
                                  ncore = 1,
                                  log_enrichment = TRUE){

  null_method <- match.arg(null_method)

  # estimated enrichment
  enrichment <- group_prior[names(group_prior) != "SNP"]/group_prior[names(group_prior) == "SNP"]

  loglik <- compute_loglik_ser(group_prior = group_prior,
                               group_prior_var = group_prior_var,
                               region_data = region_data,
                               null_method = null_method)

  # null model (enrichment = 1)
  EM_null_res <- fit_EM(region_data,
                        niter = niter,
                        init_group_prior = group_prior,
                        init_group_prior_var = group_prior_var,
                        group_prior_var_structure = "fixed",
                        shared_group_prior = TRUE,
                        null_method = null_method,
                        ncore = ncore)
  group_prior_null <- EM_null_res$group_prior
  group_prior_var_null <- EM_null_res$group_prior_var
  loglik_null <- EM_null_res$loglik_iters[length(EM_null_res$loglik_iters)]

  enrichment_null <- group_prior_null[names(group_prior_null) != "SNP"]/group_prior_null[names(group_prior_null) == "SNP"]
  stopifnot(all(enrichment_null == 1))

  if (log_enrichment){
    enrichment <- log(enrichment)
    enrichment_null <- 0
  }

  # compute S.E. using asymptotic test based on LRT
  enrichment.se <- sqrt( enrichment^2 / (2 * (loglik - loglik_null)) )

  # p-value of enrichment
  enrichment.z <- (enrichment - enrichment_null) / enrichment.se
  enrichment.pval <- 2*pnorm(abs(enrichment.z), lower.tail=FALSE)

  return(list("enrichment" = enrichment,
              "se" = enrichment.se,
              "p.value" = enrichment.pval))
}

#' Compute standard error and p-value for enrichment using the G test
#'
#' @param region_data a list of assembled region data.
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#'
#' @param group_prior_var a vector of prior variances for different groups.
#'
#' @param null_method Method to compute null model, options: "ctwas", "susie" or "none".
#'
#' @param ncore The number of cores used to parallelize over regions
#'
#' @return Estimated enrichment, S.E. and p-value.
#'
#' @export
#'
#' @importFrom parallel mclapply
#' @importFrom AMR g.test
#'
get_enrichment_se_G_test <- function(region_data,
                                     group_prior = NULL,
                                     group_prior_var = NULL,
                                     null_method = c("ctwas", "susie", "none"),
                                     ncore = 1){

  # get groups, types and contexts from region_data
  groups <- unique(unlist(lapply(region_data, "[[", "groups")))
  groups <- c(setdiff(groups, "SNP"), "SNP")

  # set pi_prior and V_prior based on group_prior and group_prior_var
  res <- initiate_group_priors(group_prior[groups], group_prior_var[groups], groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  # run finemapping with SER model to get PIPs
  region_ids <- names(region_data)
  all_ser_res_list <- mclapply_check(region_ids, function(region_id){
    finemap_single_region_ser_rss(region_data, region_id, pi_prior, V_prior,
                                  null_method = null_method,
                                  return_full_result = TRUE)
  }, mc.cores = ncore, stop_if_missing = TRUE)
  ser_res_df <- do.call(rbind, lapply(all_ser_res_list, "[[", "ser_res_df"))

  group_pip <- sapply(groups, function(x){sum(ser_res_df$susie_pip[ser_res_df$group==x])})
  names(group_pip) <- groups

  group_size <- table(ser_res_df$group)
  group_size <- group_size[groups]
  group_size <- as.numeric(group_size)
  names(group_size) <- groups

  gene_groups <- setdiff(groups, "SNP")

  enrichment <- rep(NA, length(gene_groups))
  names(enrichment) <- gene_groups
  enrichment.se <- rep(NA, length(gene_groups))
  names(enrichment.se) <- gene_groups
  enrichment.pval <- rep(NA, length(gene_groups))
  names(enrichment.pval) <- gene_groups

  for(group in gene_groups){
    k0 <- group_size["SNP"]
    k1 <- group_size[group]

    r0 <- group_pip["SNP"]
    r1 <- group_pip[group]

    # enrichment: relative risk estimate in a 2 × 2 contingency table
    enrichment[group] <- log((r1 * k0) / (r0 * k1))

    # standard error: standard error of a relative risk
    enrichment.se[group] <- sqrt(1/r1 + 1/r0 - 1/k1 - 1/k0)

    obs <- rbind(c(r1, k1),
                 c(r0, k0))
    colnames(obs) <- c("r", "k")
    rownames(obs) <- c(group, "SNP")

    # evaluate statistical significance of enrichment using the G-test
    G_test_res <- suppressWarnings(g.test(obs))

    enrichment.pval[group] <- G_test_res$p.value
  }

  return(list("enrichment" = enrichment,
              "se" = enrichment.se,
              "p.value" = enrichment.pval))
}
