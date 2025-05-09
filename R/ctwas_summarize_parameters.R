#' @title Summarizes estimated parameters
#'
#' @param param a list of parameter estimation result from \code{est_param}
#'
#' @param gwas_n the sample size of the GWAS summary statistics.
#'
#' @param enrichment_test Method to test enrichment,
#'"G": G-test, "fisher": Fisher's exact test,
#'"z": p-value computed from z-score of enrichment.
#'
#' @param alternative indicates the alternative hypothesis and
#' must be one of "greater", "two.sided", or "less".
#' Only used when \code{enrichment_test} is "fisher" or "z".
#'
#' @param include_test_result If TRUE, return the original test result.
#'
#' @return a list of summarized parameters
#'
#' @export
#'
summarize_param <- function(param,
                            gwas_n,
                            enrichment_test = c("fisher","G", "z"),
                            alternative = c("greater","two.sided","less"),
                            include_test_result = FALSE){

  enrichment_test <- match.arg(enrichment_test)
  alternative <- match.arg(alternative)

  # estimated group prior inclusion probabilieis
  group_prior <- param$group_prior

  # estimated group prior variance
  group_prior_var <- param$group_prior_var

  # set group size
  group_size <- param$group_size
  group_size <- group_size[names(group_prior)]

  # compute enrichment, s.e. and p-values
  enrichment_res <- compute_enrichment_test(group_prior = group_prior,
                                            group_size = group_size,
                                            enrichment_test = enrichment_test,
                                            alternative = alternative,
                                            include_test_result = include_test_result)
  enrichment <- enrichment_res$enrichment
  enrichment_se <- enrichment_res$se
  enrichment_pval <- enrichment_res$p.value
  enrichment_test_res <- enrichment_res$test_res

  # estimated group PVE
  group_pve <- group_prior_var*group_prior*group_size/gwas_n

  # total PVE
  total_pve <- sum(group_pve, na.rm = TRUE)

  # proportion of heritability
  prop_heritability <- group_pve/total_pve

  return(list("group_size" = group_size,
              "group_prior" = group_prior,
              "group_prior_var" = group_prior_var,
              "enrichment" = enrichment,
              "enrichment_se" = enrichment_se,
              "enrichment_pval" = enrichment_pval,
              "enrichment_test_res" = enrichment_test_res,
              "group_pve" = group_pve,
              "total_pve" = total_pve,
              "prop_heritability" = prop_heritability))
}


#' @title Computes enrichment (log-scale), standard error and p-value
#'
#' @param group_prior a vector of prior inclusion probabilities for different groups.
#'
#' @param group_size a vector of number of variables in different groups.
#'
#' @param enrichment_test Method to test enrichment,
#'"G": G-test, "fisher": Fisher's exact test,
#'"z": p-value computed from z-score of enrichment.
#'
#' @param alternative indicates the alternative hypothesis and
#' must be one of "greater", "two.sided", or "less".
#' Only used when \code{enrichment_test} is "fisher" or "z".
#'
#' @param include_test_result If TRUE, return the original test result.
#'
#' @return Estimated enrichment, S.E. and p-value from G-test or Fisher's exact test.
#'
#' @importFrom AMR g.test
#' @importFrom stats fisher.test
#'
#' @export
#'
compute_enrichment_test <- function(group_prior,
                                    group_size,
                                    enrichment_test = c("fisher","G", "z"),
                                    alternative = c("greater","two.sided","less"),
                                    include_test_result = FALSE){

  enrichment_test <- match.arg(enrichment_test)
  alternative <- match.arg(alternative)

  # compute sum of PIPs for each group
  groups <- names(group_prior)
  group_size <- group_size[groups]
  group_sum_pip <- group_prior * group_size
  names(group_sum_pip) <- groups

  gene_groups <- setdiff(groups, "SNP")

  enrichment <- rep(NA, length(gene_groups))
  names(enrichment) <- gene_groups
  enrichment.se <- rep(NA, length(gene_groups))
  names(enrichment.se) <- gene_groups
  enrichment.pval <- rep(NA, length(gene_groups))
  names(enrichment.pval) <- gene_groups

  test_res_list <- list()
  for (group in gene_groups) {
    k0 <- group_size["SNP"]
    k1 <- group_size[group]

    r0 <- group_sum_pip["SNP"]
    r1 <- group_sum_pip[group]

    # enrichment: relative risk estimate in a 2 × 2 contingency table
    enrichment[group] <- log((r1/k1) / (r0/k0))

    # standard error: standard error of a relative risk
    enrichment.se[group] <- sqrt(1/r1 + 1/r0 - 1/k1 - 1/k0)

    obs <- rbind(c(r1, k1), c(r0, k0))
    colnames(obs) <- c("r", "k")
    rownames(obs) <- c(group, "SNP")

    if (enrichment_test == "G"){
      # evaluate statistical significance of enrichment using G-test
      test_res <- g.test(obs)
      enrichment.pval[group] <- test_res$p.value
    } else if (enrichment_test == "fisher"){
      # evaluate statistical significance of enrichment using Fisher's exact test
      test_res <- fisher.test(round(obs), alternative = alternative)
      enrichment.pval[group] <- test_res$p.value
    } else if (enrichment_test == "z"){
      enrichment <- as.numeric(enrichment[group])
      se <- as.numeric(enrichment.se[group])
      z <- enrichment/se
      if (alternative == "two.sided"){
        p <- 2*pnorm(abs(z), lower.tail=FALSE)
      } else if (alternative == "greater") {
        p <- pnorm(abs(z), lower.tail=FALSE)
      } else if (alternative == "less") {
        p <- pnorm(-abs(z), lower.tail=TRUE)
      }
      enrichment.pval[group] <- p
      test_res <- list("enrichment" = enrichment, "se" = se, "z" = z, "p" = p)
    }

    # save test result
    test_res_list[[group]] <- test_res
  }

  if (!include_test_result) {
    test_res_list <- NULL
  }

  return(list("enrichment" = enrichment,
              "se" = enrichment.se,
              "p.value" = enrichment.pval,
              "test_res" = test_res_list))
}
