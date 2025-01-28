# fit cTWAS version of single effect regression (SER) model with summary statistics
ctwas_ser_rss <- function(z,
                          prior_variance,
                          prior_weights = NULL,
                          residual_variance = 1,
                          null_method = c("ctwas", "susie", "none"),
                          normalize_prior = TRUE){

  null_method = match.arg(null_method)

  if (length(prior_weights) != length(z))
    stop("Prior weights must have length of z")

  if (any(prior_weights < 0 | prior_weights > 1))
    warning("Prior weights must be between 0 and 1")

  p = length(z)
  shat2 = 1
  betahat = z

  if (is.null(prior_weights))
    prior_weights = rep(1/p,p)

  if (null_method == "susie") {
    null_weight = max(0, 1 - sum(prior_weights))
    if (normalize_prior){
      # rescale prior_weights and null_weight, so that they sum to 1
      sum_prior_weights = sum(prior_weights) + null_weight
      prior_weights = prior_weights / sum_prior_weights
      null_weight = null_weight / sum_prior_weights
    }
    prior_null = null_weight
    prior_configs = prior_weights
  } else if (null_method == "ctwas") {
    prior_null = prod(1 - prior_weights)
    prior_odds = prior_weights / (1-prior_weights)
    prior_configs = prior_odds * prior_null
  } else if (null_method == "none") {
    prior_null = 0
    if (normalize_prior){
      # rescale prior_weights, so that they sum to 1
      prior_weights = prior_weights / sum(prior_weights)
    }
    prior_configs = prior_weights
  }

  # compute BF using the Wakefield's formula
  gamma = prior_variance / (prior_variance + shat2)
  lbf = log(sqrt(1-gamma)) - (-1/2 * z^2  * gamma)
  maxlbf = max(lbf)
  maxbf = exp(maxlbf)

  # w is proportional to BF, but subtract max for numerical stability.
  w = exp(lbf - maxlbf)

  # posterior prob for each variable
  null_w = prior_null/maxbf
  weighted_w = prior_configs * w
  weighted_sum_w = sum(weighted_w)

  if (null_method == "ctwas") {
    weighted_w_scaled = prior_odds * w
    alpha = weighted_w_scaled / (1/maxbf + sum(weighted_w_scaled))
  } else {
    alpha = weighted_w / (null_w + weighted_sum_w)
  }

  # compute posterior mean and second moment
  post_var = (1/prior_variance + 1/shat2)^(-1) # posterior variance
  post_mean =  post_var * betahat / shat2 # posterior mean
  post_mean2 = post_var + post_mean^2 # Second moment

  # compute the marginal log-likelihood (minus the constant of loglik_null)
  loglik = log(null_w + weighted_sum_w) + maxlbf

  return(list(alpha = alpha,
              post_var = post_var,
              mu = post_mean,
              mu2 = post_mean2,
              lbf = lbf,
              loglik = loglik,
              prior_null = prior_null,
              prior_configs = prior_configs,
              prior_variance = prior_variance))
}

# fit cTWAS version of single effect regression (SER) model with summary statistics
ctwas_ser <- function(X, Y,
                      scaled_prior_variance,
                      prior_weights = NULL,
                      residual_variance = 1,
                      null_method = c("ctwas", "susie", "none"),
                      normalize_prior = TRUE){

  null_method = match.arg(null_method)

  if (length(prior_weights) != ncol(X))
    stop("Prior weights must have length of ncol(X)")

  if (any(prior_weights < 0 | prior_weights > 1))
    warning("Prior weights must be between 0 and 1")

  p = ncol(X)
  n = nrow(X)

  # Center and scale input.
  X = set_X_attributes(X, center = TRUE, scale = TRUE)
  Y = scale(Y, center = TRUE, scale = TRUE)

  XtX = crossprod(X)
  Xty = drop(crossprod(X,y))
  dXtX = diag(XtX)

  betahat = (1/dXtX) * Xty
  shat2 = residual_variance / dXtX
  z = betahat / sqrt(shat2)

  if (is.null(prior_weights))
    prior_weights = rep(1/p,p)

  if (null_method == "susie") {
    null_weight = max(0, 1 - sum(prior_weights))
    if (normalize_prior){
      # rescale prior_weights and null_weight, so that they sum to 1
      sum_prior_weights = sum(prior_weights) + null_weight
      prior_weights = prior_weights / sum_prior_weights
      null_weight = null_weight / sum_prior_weights
    }
    prior_null = null_weight
    prior_configs = prior_weights
  } else if (null_method == "ctwas") {
    prior_null = prod(1 - prior_weights)
    prior_odds = prior_weights / (1-prior_weights)
    prior_configs = prior_odds * prior_null
  } else if (null_method == "none") {
    prior_null = 0
    if (normalize_prior){
      # rescale prior_weights, so that they sum to 1
      prior_weights = prior_weights / sum(prior_weights)
    }
    prior_configs = prior_weights
  }

  # compute BF using the Wakefield's formula
  gamma = scaled_prior_variance / (scaled_prior_variance + shat2)
  lbf = log(sqrt(1-gamma)) - (-1/2 * z^2  * gamma)
  maxlbf = max(lbf)
  maxbf = exp(maxlbf)

  # w is proportional to BF, but subtract max for numerical stability.
  w = exp(lbf - maxlbf)

  # posterior prob for each variable
  null_w = prior_null/maxbf
  weighted_w = prior_configs * w
  weighted_sum_w = sum(weighted_w)

  if (null_method == "ctwas") {
    weighted_w_scaled = prior_odds * w
    alpha = weighted_w_scaled / (1/maxbf + sum(weighted_w_scaled))
  } else {
    alpha = weighted_w / (null_w + weighted_sum_w)
  }

  # compute posterior mean and second moment
  post_var = (1/scaled_prior_variance + 1/shat2)^(-1) # posterior variance
  post_mean =  post_var * betahat / shat2 # posterior mean
  post_mean2 = post_var + post_mean^2 # Second moment

  # compute the marginal log-likelihood (minus the constant of loglik_null)
  loglik = log(null_w + weighted_sum_w) + maxlbf

  return(list(alpha = alpha,
              post_var = post_var,
              mu = post_mean,
              mu2 = post_mean2,
              lbf = lbf,
              loglik = loglik,
              prior_null = prior_null,
              prior_configs = prior_configs,
              scaled_prior_variance = scaled_prior_variance))
}


# annotate SER results with gene and SNP information
anno_ser_res <- function(ser_res,
                         gids,
                         sids,
                         g_type = "gene",
                         g_context = "gene",
                         g_group = "gene",
                         region_id = NULL,
                         z = NULL) {

  gene_df <- data.frame(id = gids,
                        type = g_type, context = g_context, group = g_group)

  snp_df <- data.frame(id = sids,
                       type = "SNP", context = "SNP", group = "SNP")

  ser_res_df <- rbind(gene_df, snp_df)

  if (!is.null(region_id)) {
    ser_res_df$region_id <- region_id
  }

  if (!is.null(z)) {
    ser_res_df$z <- z
  }

  ser_res_df$susie_pip <- drop(ser_res$alpha)
  ser_res_df$mu2 <- drop(ser_res$mu2)

  return(ser_res_df)

}

# finemap a single region with L = 1 without LD, used in EM
finemap_single_region_ser_rss <- function(region_data,
                                          region_id,
                                          pi_prior,
                                          V_prior,
                                          null_method = c("ctwas", "susie", "none"),
                                          return_full_result = FALSE){

  null_method <- match.arg(null_method)

  # load region data
  if (!inherits(region_data,"list")){
    stop("'region_data' should be a list.")
  }

  regiondata <- region_data[[region_id]]
  if (is.null(regiondata$z) || is.null(regiondata$gs_group)){
    regiondata <- extract_region_data(region_data, region_id)
  }

  gids <- regiondata[["gid"]]
  sids <- regiondata[["sid"]]
  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  g_type <- regiondata[["g_type"]]
  g_context <- regiondata[["g_context"]]
  g_group <- regiondata[["g_group"]]
  rm(regiondata)

  if (length(z) < 2) {
    stop(paste(length(z), "variables in the region", region_id, "\n",
               "At least two variables in a region are needed to run SER model"))
  }

  # update priors, prior variances and null_weight
  res <- set_region_susie_priors(pi_prior, V_prior, gs_group, L = 1,
                                 use_null_weight = FALSE)
  prior_weights <- res$prior
  prior_variance <- res$V
  rm(res)

  # fit SER model
  ser_res <- ctwas_ser_rss(z = z,
                           prior_weights = prior_weights,
                           prior_variance = prior_variance,
                           null_method = null_method)

  # annotate SER result
  ser_res_df <- anno_ser_res(ser_res,
                             gids,
                             sids,
                             g_type = g_type,
                             g_context = g_context,
                             g_group = g_group,
                             region_id = region_id,
                             z = z)

  if (return_full_result){
    return(list("ser_res" = ser_res,
                "ser_res_df" = ser_res_df))
  }else{
    return(ser_res_df)
  }

}

# fit a single region with L = 1 without LD, used in EM
fit_single_region_ser_rss <- function(region_data,
                                      region_id,
                                      pi_prior,
                                      V_prior,
                                      null_method = c("ctwas", "susie", "none")){

  null_weight_method <- match.arg(null_weight_method)

  # load region data
  if (!inherits(region_data,"list")){
    stop("'region_data' should be a list.")
  }

  regiondata <- region_data[[region_id]]
  if (is.null(regiondata$z) || is.null(regiondata$gs_group)){
    regiondata <- extract_region_data(region_data, region_id)
  }

  z <- regiondata[["z"]]
  gs_group <- regiondata[["gs_group"]]
  rm(regiondata)

  if (length(z) < 2) {
    stop(paste(length(z), "variables in the region", region_id, "\n",
               "At least two variables in a region are needed to run SER model"))
  }

  # update priors, prior variances and null_weight
  res <- set_region_susie_priors(pi_prior, V_prior, gs_group, L = 1,
                                 use_null_weight = FALSE)
  prior_weights <- res$prior
  prior_variance <- res$V
  rm(res)

  # fit SER model
  ser_res <- ctwas_ser_rss(z = z,
                           prior_weights = prior_weights,
                           prior_variance = prior_variance,
                           null_method = null_method)

  return(ser_res)
}


#' @importFrom logging addHandler loginfo logwarn writeToFile
#'
compute_loglik_ser <- function(
    group_prior,
    group_prior_var,
    region_data,
    groups,
    null_method = c("ctwas", "susie", "none"),
    ncore = 1,
    verbose = FALSE){

  # get groups, from region_data
  if (missing(groups)){
    groups <- unique(unlist(lapply(region_data, "[[", "groups")))
    groups <- c(setdiff(groups, "SNP"), "SNP")
  }

  if (verbose) {
    loginfo("group_prior {%s}: {%s}",
            names(group_prior), format(group_prior, digits = 4))
    loginfo("group_prior_var {%s}: {%s}",
            names(group_prior_var), format(group_prior_var, digits = 4))
  }

  res <- initiate_group_priors(group_prior[groups], group_prior_var[groups], groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  region_ids <- names(region_data)

  EM_ser_res_list <- mclapply_check(region_ids, function(region_id){
    finemap_single_region_ser_rss(region_data, region_id, pi_prior, V_prior,
                                  null_method = null_method,
                                  return_full_result = TRUE)
  }, mc.cores = ncore, stop_if_missing = TRUE)

  all_ser_res <- lapply(EM_ser_res_list, "[[", "ser_res")
  loglik <- sum(sapply(all_ser_res, "[[", "loglik"))

  return(loglik)
}
