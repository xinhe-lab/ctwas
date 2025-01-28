
# Run cTWAS version of susie_rss for a single region
ctwas_susie_rss <- function(z,
                            R,
                            prior_weights = NULL,
                            prior_variance = NULL,
                            L = 5,
                            z_ld_weight = 0,
                            null_weight = NULL,
                            coverage = 0.95,
                            min_abs_corr = 0.1,
                            warn_converge_fail = TRUE,
                            ...){

  if (missing(R)) {
    if (L == 1) {
      # R does not matter for susie when L = 1
      R <- diag(length(z))
    } else {
      stop("R (correlation matrix) is required when L > 1")
    }
  }

  if (L == 1) {
    # turn off purity cutoff when L = 1
    min_abs_corr <- 0
  }

  # in susie, prior_variance is under standardized scale (if performed)
  susie_res <- susie_rss(z,
                         R,
                         prior_weights = prior_weights,
                         prior_variance = prior_variance,
                         estimate_prior_variance = FALSE,
                         L = L,
                         z_ld_weight = z_ld_weight,
                         null_weight = null_weight,
                         coverage = coverage,
                         min_abs_corr = min_abs_corr,
                         warn_converge_fail = warn_converge_fail,
                         ...)

  return(susie_res)
}


# annotate susie results with SNP and gene information
anno_susie <- function(susie_res,
                       gids,
                       sids,
                       g_type = "gene",
                       g_context = "gene",
                       g_group = "gene",
                       region_id = NULL,
                       z = NULL,
                       include_cs = TRUE,
                       get_alpha = TRUE) {

  gene_df <- data.frame(id = gids,
                        molecular_id = sapply(strsplit(gids, split = "[|]"), "[[", 1),
                        type = g_type, context = g_context, group = g_group)

  snp_df <- data.frame(id = sids,
                       molecular_id = sids,
                       type = "SNP", context = "SNP", group = "SNP")

  susie_res_df <- rbind(gene_df, snp_df)

  if (!is.null(region_id)) {
    susie_res_df$region_id <- region_id
  }

  if (!is.null(z)) {
    susie_res_df$z <- z
  }

  susie_res_df$susie_pip <- susie_res$pip

  p <- length(gids) + length(sids)
  susie_res_df$mu2 <- colSums(susie_res$mu2[, seq(1, p)[1:p!=susie_res$null_index], drop = FALSE]) # WARN: not sure for L>1

  if (include_cs) {
    susie_res_df$cs <- NA
    if (!is.null(susie_res$sets$cs)){
      for (cs_i in susie_res$sets$cs_index){
        cs_name <- paste0("L", cs_i)
        tmp.cs.idx <- susie_res$sets$cs[[cs_name]]

        # Drop null weight indices
        if (!is.null(susie_res$null_index) && susie_res$null_index > 0){
          tmp.cs.idx <- tmp.cs.idx[tmp.cs.idx != susie_res$null_index]
        }

        # append cs
        susie_res_df$cs[tmp.cs.idx] <- ifelse(
          is.na(susie_res_df$cs[tmp.cs.idx]),
          cs_name,
          paste0(susie_res_df$cs[tmp.cs.idx], ",", cs_name))
      }
    }
  }

  return(susie_res_df)

}

#' extract the alpha matrix from susie result
#' adapted from susie's susie_get_pip() function
#'
#' @param prune_by_cs Whether or not to ignore single effects not in
#'   a reported CS when calculating PIP.
#'
#' @keywords internal
#'
extract_susie_alpha = function (susie_res, prune_by_cs = FALSE) {

  if (inherits(susie_res,"susie")) {

    # Drop null weight columns.
    if (!is.null(susie_res$null_index) && susie_res$null_index > 0)
      susie_res$alpha = susie_res$alpha[,-susie_res$null_index,drop=FALSE]

    include_idx = 1:nrow(susie_res$alpha)

    # Only consider variables in reported CS.
    # This is not what we do in the SuSiE paper.
    # So by default prune_by_cs = FALSE means we do not run the
    # following code.
    if (!is.null(susie_res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx,susie_res$sets$cs_index)
    if (is.null(susie_res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)

    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0)
      susie_alpha = susie_res$alpha[include_idx,,drop = FALSE]
    else
      susie_alpha = matrix(0,1,ncol(susie_res$alpha))
  }

  return(susie_alpha)
}

# extract the alpha matrix from susie result,
# and combine with annotated susie result data frame
#' @importFrom tidyr pivot_longer
get_susie_alpha_res <- function(susie_res,
                                susie_res_df,
                                keep_genes_only = TRUE) {

  # extract alpha matrix from susie result
  susie_alpha <- as.data.frame(t(extract_susie_alpha(susie_res)))
  colnames(susie_alpha) <- paste0("L", seq_len(ncol(susie_alpha)))
  rownames(susie_alpha) <- susie_res_df$id

  susie_alpha_df <- cbind(susie_res_df, susie_alpha)

  susie_alpha_df <- pivot_longer(susie_alpha_df,
                                 cols = colnames(susie_alpha),
                                 names_to = "susie_set",
                                 values_to = "susie_alpha")
  susie_alpha_df <- as.data.frame(susie_alpha_df)

  susie_alpha_df$in_cs <- FALSE
  if (!is.null(susie_res$sets$cs)){
    for (cs_i in susie_res$sets$cs_index){
      cs_name <- paste0("L", cs_i)
      tmp.cs.idx <- susie_res$sets$cs[[cs_name]]

      # Drop null weight indices
      if (!is.null(susie_res$null_index) && susie_res$null_index > 0){
        tmp.cs.idx <- tmp.cs.idx[tmp.cs.idx != susie_res$null_index]
      }

      tmp.cs.id <- susie_res_df$id[tmp.cs.idx]

      # flag the alpha in cs
      susie_alpha_df$in_cs[which(susie_alpha_df$susie_set == cs_name & susie_alpha_df$id %in% tmp.cs.id)] <- TRUE
    }
  }

  if (keep_genes_only) {
    susie_alpha_df <- susie_alpha_df[susie_alpha_df$group != "SNP",,drop=FALSE]
  }

  return(susie_alpha_df)
}

# set pi_prior and V_prior based on init_group_prior and init_group_prior_var
initiate_group_priors <- function(group_prior = NULL, group_prior_var = NULL, groups) {

  if (is.null(group_prior)){
    group_prior <- structure(as.numeric(rep(NA,length(groups))), names=groups)
  }

  if (is.null(group_prior_var)){
    group_prior_var <- structure(as.numeric(rep(NA,length(groups))), names=groups)
  }

  if ( !setequal(names(group_prior), groups) || !setequal(names(group_prior_var), groups) ) {
    stop("Names of group_prior or group_prior_var do not match with groups")
  }

  pi_prior <- list()
  V_prior <- list()
  for (group in groups){
    pi_prior[[group]] <- unname(group_prior[group])
    V_prior[[group]] <- unname(group_prior_var[group])
  }
  pi_prior <- unlist(pi_prior)
  V_prior <- unlist(V_prior)

  return(list("pi_prior" = pi_prior,
              "V_prior" = V_prior))
}


# set prior and prior variance values for the region
set_region_susie_priors <- function(pi_prior, V_prior, gs_group, L,
                                    use_null_weight = TRUE){

  if (length(gs_group) < 2) {
    stop(paste(length(gs_group), "variables in the region. At least two variables in a region are needed to run susie"))
  }

  p <- length(gs_group)

  if (any(is.na(pi_prior))){
    prior <- rep(1/p, p)
  } else {
    prior <- unname(pi_prior[gs_group])
  }

  if (any(is.na(V_prior))){
    # following the default in susie_rss of susieR
    V <- matrix(rep(50, L * p), nrow = L)
  } else{
    V <- unname(V_prior[gs_group])
    V <- matrix(rep(V, each = L), nrow=L)
  }

  if (use_null_weight){
    null_weight <- max(0, 1 - sum(prior))
    prior <- prior/(1-null_weight)
  } else {
    null_weight <- NULL
  }

  return(list("prior" = prior,
              "V" = V,
              "null_weight" = null_weight))

}
