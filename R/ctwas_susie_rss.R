
# Run cTWAS version of susie_rss for a single region
ctwas_susie_rss <- function(z,
                            R,
                            prior_weights = NULL,
                            prior_variance = NULL,
                            L = 5,
                            z_ld_weight = 0,
                            null_weight = NULL,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            warn_converge_fail = TRUE,
                            ...){

  if (missing(R)) {
    if (L == 1){
      # R does not matter for susie when L = 1
      R <- diag(length(z))
    } else {
      stop("R (correlation matrix) is required when L > 1")
    }
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
                       region_id,
                       z = NULL,
                       g_type = "gene",
                       g_context = "gene",
                       g_group = "gene",
                       include_cs_index = TRUE) {

  gene_df <- data.frame(id = gids, type = g_type, context = g_context, group = g_group)
  snp_df <- data.frame(id = sids, type = "SNP", context = "SNP", group = "SNP")
  susie_res_df <- rbind(gene_df, snp_df)

  if (!is.null(z)) {
    susie_res_df$z <- z
  }

  susie_res_df$region_id <- region_id

  susie_res_df$susie_pip <- susie_res$pip

  p <- length(gids) + length(sids)
  susie_res_df$mu2 <- colSums(susie_res$mu2[, seq(1, p)[1:p!=susie_res$null_index], drop = F]) #WARN: not sure for L>1

  if (include_cs_index) {
    susie_res_df$cs_index <- 0
    if (!is.null(susie_res$sets$cs)){
      for (cs_i in susie_res$sets$cs_index){
        X.idx <- susie_res$sets$cs[[paste0("L", cs_i)]]
        X.idx <- X.idx[X.idx != susie_res$null_index] # susie_rss' bug
        susie_res_df$cs_index[X.idx] <- cs_i
        #TODO: note this ignores the fact that some variants can belong to multiple CS
      }
    }
  }

  return(susie_res_df)
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
    stop("names of group_prior or group_prior_var do not match with groups")
  }

  pi_prior <- list()
  V_prior <- list()
  for (group in groups){
    pi_prior[[group]] <- unname(group_prior[group])
    V_prior[[group]] <- unname(group_prior_var[group])
  }
  pi_prior <- unlist(pi_prior)
  V_prior <- unlist(V_prior)

  return(list(pi_prior = pi_prior,
              V_prior = V_prior))
}


# set prior and prior variance values for the region
set_region_susie_priors <- function(pi_prior, V_prior, gs_group, L, use_null_weight = TRUE){

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
    V <- matrix(rep(50, L * p), nrow = L)
    # following the default in susie_rss of susieR
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

  return(list(prior = prior,
              V = V,
              null_weight = null_weight))

}
