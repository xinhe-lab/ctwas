# simulate a single group
simulate_single_group <- function(n,
                                  p,
                                  prior_weights,
                                  prior_variance,
                                  effect_scale = c("z", "beta"),
                                  simulate_null_model = FALSE,
                                  max_n_effect = Inf,
                                  min_n_effect = 0,
                                  compute_R = FALSE,
                                  keep_sumstats_only = FALSE,
                                  seed = 1,
                                  verbose = FALSE){

  effect_scale <- match.arg(effect_scale)

  set.seed(seed)

  if (effect_scale == "z")
    prior_variance <- prior_variance / n

  # simulate beta
  config_res <- simulate_config_single_group(p, prior_weights, simulate_null_model)
  # resampled more then one effect, resimulate
  while(length(config_res$effect.idx) > max_n_effect | length(config_res$effect.idx) < min_n_effect){
    if(verbose)
      cat("Sampled", length(config_res$effect.idx), "causal effects, resampling...\n")
    config_res <- simulate_config_single_group(p, prior_weights, simulate_null_model)
  }
  if(verbose)
    cat("Sampled", length(config_res$effect.idx), "causal effects. \n")
  gamma <- config_res$gamma
  effect.idx <- config_res$effect.idx

  beta = rep(0,p)
  if (length(effect.idx) > 0) {
    beta[effect.idx] = rnorm(length(effect.idx), mean = 0, sd = sqrt(prior_variance[effect.idx]))
  }
  # simulate X, y, z
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  X = scale(X, center = TRUE, scale = TRUE)
  y = drop(X %*% beta + rnorm(n))
  y = drop(scale(y,center = TRUE,scale = TRUE))

  XtX = crossprod(X)
  Xty = drop(crossprod(X,y))
  yty = sum(y^2)
  dXtX = diag(XtX)
  betahat = (1/dXtX) * Xty
  shat = sqrt(1/dXtX)
  z = betahat / shat

  # input_ss = compute_ss(X,y,standardize = TRUE)
  # R    = with(input_ss,cov2cor(XtX))
  # ss = univariate_regression(X,y)
  # z = ss$betahat/ss$sebetahat

  if (compute_R){
    R = cov2cor(XtX)
  } else {
    R = NULL
  }

  if (keep_sumstats_only) {
    X <- NULL
    y <- NULL
    beta <- NULL
    gamma <- NULL
  }

  return(list(n = n,
              p = p,
              X = X,
              y = y,
              z = z,
              R = R,
              beta = beta,
              gamma = gamma,
              effect.idx = effect.idx))
}

# simulate two groups (genes and snps)
simulate_groups <- function(group_size,
                            n,
                            group_prior = NULL,
                            group_prior_var = NULL,
                            effect_scale = c("z", "beta"),
                            simulate_null_model = FALSE,
                            max_n_effect = Inf,
                            min_n_effect = 0,
                            compute_R = FALSE,
                            keep_sumstats_only = FALSE,
                            seed = 1,
                            verbose = FALSE){

  effect_scale <- match.arg(effect_scale)

  set.seed(seed)
  p = sum(group_size)
  groups = names(group_size)

  # convert to beta scale
  if (effect_scale == "z")
    group_prior_var <- group_prior_var / n

  res <- initialize_group_priors(group_prior[groups], group_prior_var[groups], groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  if (max_n_effect == 0)
    simulate_null_model = TRUE

  config_res <- simulate_config_groups(group_size, pi_prior, simulate_null_model)
  # resampled more then one effect, resimulate
  while(length(config_res$effect.idx) > max_n_effect | length(config_res$effect.idx) < min_n_effect){
    if(verbose)
      cat("Sampled", length(config_res$effect.idx), "causal effects, resampling...\n")
    config_res <- simulate_config_groups(group_size, pi_prior, simulate_null_model)
  }
  if(verbose)
    cat("Sampled", length(config_res$effect.idx), "causal effects. \n")
  gamma <- config_res$gamma
  gidx <- config_res$gidx
  sidx <- config_res$sidx
  effect.idx <- config_res$effect.idx
  gid <- config_res$gid
  sid <- config_res$sid
  gs_group <- config_res$gs_group

  # set prior and prior variance values for the variables
  if (any(is.na(pi_prior))){
    prior_weights <- rep(1/p, p)
  } else {
    prior_weights <- unname(pi_prior[gs_group])
  }

  if (any(is.na(V_prior))){
    prior_variance_z <- matrix(rep(50/n, p), nrow = 1)
    prior_variance_beta <- prior_variance_z / n
  } else{
    prior_variance_beta <- matrix(unname(V_prior[gs_group]), nrow=1)
    prior_variance_z = prior_variance_beta * n
  }

  # simulate beta
  beta = rep(0,p)
  if (length(effect.idx) > 0) {
    beta[effect.idx] = rnorm(length(effect.idx), mean = 0, sd = sqrt(drop(prior_variance_beta)[effect.idx]))
  }
  # simulate X, y, z
  X = matrix(rnorm(n*p),nrow = n,ncol = p)
  X = scale(X,center = TRUE,scale = TRUE)
  y = drop(X %*% beta + rnorm(n))
  y = drop(scale(y,center = TRUE,scale = TRUE))

  XtX = crossprod(X)
  Xty = drop(crossprod(X,y))
  yty = sum(y^2)
  dXtX = diag(XtX)
  betahat = (1/dXtX) * Xty
  shat = sqrt(1/dXtX)
  z = betahat / shat

  # input_ss = compute_ss(X,y,standardize = TRUE)
  # R    = with(input_ss,cov2cor(XtX))
  # ss = univariate_regression(X,y)
  # zhat = ss$betahat/ss$sebetahat

  if (compute_R){
    R = cov2cor(XtX)
  } else {
    R = NULL
  }

  if (keep_sumstats_only) {
    X <- NULL
    y <- NULL
    beta <- NULL
    gamma <- NULL
  }

  return(list(n = n,
              p = p,
              X = X,
              y = y,
              z = z,
              R = R,
              beta = beta,
              gamma = gamma,
              effect.idx = effect.idx,
              gs_group = gs_group,
              gid = gid,
              sid = sid,
              group_size = group_size,
              groups = groups,
              prior_weights = prior_weights,
              prior_variance_beta = prior_variance_beta,
              prior_variance_z = prior_variance_z))
}

# simulate region data
simulate_region_data <- function(group_size,
                                 n,
                                 region_id,
                                 group_prior = NULL,
                                 group_prior_var = NULL,
                                 effect_scale = c("z", "beta"),
                                 groups,
                                 group_types,
                                 group_contexts,
                                 simulate_null_model = FALSE,
                                 max_n_effect = Inf,
                                 min_n_effect = 0,
                                 seed = 1,
                                 verbose = FALSE){

  effect_scale <- match.arg(effect_scale)

  set.seed(seed)
  p = sum(group_size)
  if (missing(groups)){
    groups <- c(setdiff(names(group_size), "SNP"), "SNP")
  }
  if (missing(group_types)){
    group_types <- names(group_size)
  }
  if (missing(group_contexts)){
    group_contexts <- names(group_size)
  }

  gene_groups = setdiff(names(group_size), "SNP")

  # convert to beta scale
  if (effect_scale == "z")
    group_prior_var <- group_prior_var / n

  res <- initialize_group_priors(group_prior[groups], group_prior_var[groups], groups)
  pi_prior <- res$pi_prior
  V_prior <- res$V_prior
  rm(res)

  if (max_n_effect == 0)
    simulate_null_model = TRUE

  config_res <- simulate_config_groups(group_size, pi_prior, simulate_null_model)
  # resampled more then one effect, resimulate
  while(length(config_res$effect.idx) > max_n_effect | length(config_res$effect.idx) < min_n_effect){
    if(verbose)
      cat("Sampled", length(config_res$effect.idx), "causal effect, resampling...\n")
    config_res <- simulate_config_groups(group_size, pi_prior, simulate_null_model)
  }
  if(verbose)
    cat("Sampled", length(config_res$effect.idx), "causal effect. \n")
  gamma <- config_res$gamma
  gidx <- config_res$gidx
  sidx <- config_res$sidx
  effect.idx <- config_res$effect.idx
  gid <- config_res$gid
  sid <- config_res$sid
  gs_group <- config_res$gs_group

  # set prior and prior variance values for the variables
  if (any(is.na(pi_prior))){
    prior_weights <- rep(1/p, p)
  } else {
    prior_weights <- unname(pi_prior[gs_group])
  }

  if (any(is.na(V_prior))){
    prior_variance_z <- matrix(rep(50/n, p), nrow = 1)
    prior_variance_beta <- prior_variance_z / n
  } else{
    prior_variance_beta <- matrix(unname(V_prior[gs_group]), nrow=1)
    prior_variance_z <- prior_variance_beta * n
  }

  # simulate beta
  beta = rep(0,p)
  if (length(effect.idx) > 0) {
    beta[effect.idx] = rnorm(length(effect.idx), mean = 0, sd = sqrt(drop(prior_variance_beta)[effect.idx]))
  }
  # simulate X, y, z
  X = matrix(rnorm(n*p),nrow = n,ncol = p)
  X = scale(X,center = TRUE,scale = TRUE)
  y = drop(X %*% beta + rnorm(n))
  y = drop(scale(y,center = TRUE,scale = TRUE))

  XtX = crossprod(X)
  Xty = drop(crossprod(X,y))
  yty = sum(y^2)
  dXtX = diag(XtX)
  betahat = (1/dXtX) * Xty
  shat = sqrt(1/dXtX)
  z = betahat / shat

  # input_ss = compute_ss(X,y,standardize = TRUE)
  # R    = with(input_ss,cov2cor(XtX))
  # ss = univariate_regression(X,y)
  # zhat = ss$betahat/ss$sebetahat

  z_gene <- data.frame("id" = gid, "z" = z[gidx])
  gene_group_idx <- which(names(group_size) != "SNP")
  z_gene$type <- rep(group_types[gene_group_idx], group_size[gene_group_idx])
  z_gene$context <- rep(group_contexts[gene_group_idx], group_size[gene_group_idx])
  z_gene$group <- rep(groups[gene_group_idx], group_size[gene_group_idx])

  z_snp <- data.frame("id" = sid, "z" = z[sidx],  "type" = "SNP", "context" = "SNP", "group" = "SNP")

  types <- unique(group_types)
  contexts <- unique(group_contexts)

  return(list("region_id" = region_id,
              "thin" = 1,
              "gid" = gid,
              "sid" = sid,
              "z_gene" = z_gene,
              "z_snp" = z_snp,
              "types" = types,
              "contexts" = contexts,
              "groups" = groups))
}

# simulate causal configuration
#' @importFrom stats fisher.test rbinom
simulate_config_single_group <- function(p, prior_weights, simulate_null_model) {

  # simulate beta
  gamma = rep(0, p)
  if (!simulate_null_model) {
    gamma = rbinom(p, 1, prior_weights)
    effect.idx = which(gamma == 1)
  } else {
    effect.idx = NULL
  }

  return(list("gamma" = gamma,
              "effect.idx" = effect.idx))
}


# simulate causal configuration
#' @importFrom stats fisher.test rbinom
simulate_config_groups <- function(group_size, pi_prior, simulate_null_model) {
  gene_groups = setdiff(names(group_size), "SNP")
  group_size = group_size[c(gene_groups, "SNP")]
  groups = names(group_size)
  p = sum(group_size)

  gid = unlist(lapply(1:length(gene_groups), function(i) {
    paste0(gene_groups[i], ".", seq(1, group_size[gene_groups[i]]))
  }))
  sid = paste0("SNP", seq(1, group_size["SNP"]))

  gidx = seq(1, sum(group_size[gene_groups]))
  sidx = setdiff(seq(1, p), gidx)

  gs_group = rep(groups, group_size[groups])

  if (!simulate_null_model) {
    gamma = unlist(lapply(1:length(groups), function(i) {
      rbinom(group_size[groups[i]], 1, pi_prior[groups[i]])
    }))
    effect.idx = which(gamma == 1)
  } else {
    gamma = rep(0, p)
    effect.idx = NULL
  }
  return(list("gamma" = gamma,
              "effect.idx" = effect.idx,
              "gidx" = gidx,
              "sidx" = sidx,
              "gid" = gid,
              "sid" = sid,
              "gs_group" = gs_group))
}

