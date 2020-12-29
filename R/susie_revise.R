#' @rdname update_each_effect
#' @title update each effect
#' @description Revised from original susieR::update_each_effect.
#'   Added susieR::: for all functions from susie
#' @param X an n by p matrix of regressor variables
#' @param Y an n vector of response variable
#' @param s a SuSiE fit
#' @param estimate_prior_variance boolean indicating whether to
#'   estimate prior variance
#' @param check_null_threshold float a threshold on the log scale to
#'   compare likelihood between current estimate and zero the null
update_each_effect = function (X, Y, s, estimate_prior_variance = FALSE,
                               estimate_prior_method = "optim",
                               check_null_threshold) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  L = nrow(s$alpha)
  if (L > 0)
    for (l in 1:L) {

      # Remove lth effect from fitted values.
      s$Xr = s$Xr - susieR:::compute_Xb(X,s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      R = Y - s$Xr

      res = susieR:::single_effect_regression(R,X,s$V[l,],s$sigma2,s$pi,
                                     estimate_prior_method,
                                     check_null_threshold)

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l,]      = res$V  #different from susieR
      s$lbf[l]    = res$lbf_model
      s$KL[l]     = -res$loglik +
        susieR:::SER_posterior_e_loglik(X,R,s$sigma2,res$alpha * res$mu,
                               res$alpha * res$mu2)

      s$Xr = s$Xr + susieR:::compute_Xb(X,s$alpha[l,] * s$mu[l,])
    }
  return(s)
}

assignInNamespace("update_each_effect", update_each_effect, "susieR")

#' update each effect rss (revised from susieR)
update_each_effect_rss <- function (R, z, s_init, Sigma, estimate_prior_variance = FALSE,
          estimate_prior_method = "optim", check_null_threshold = 0)
{
  if (!estimate_prior_variance)
    estimate_prior_method = "none"
  s = s_init
  L = nrow(s$alpha)
  if (L > 0)
    for (l in 1:L) {
      s$Rz = s$Rz - R %*% (s$alpha[l, ] * s$mu[l, ])
      r = z - s$Rz
      res = single_effect_regression_rss(as.vector(r),
                                         Sigma, s$V[l], s$pi, estimate_prior_method, check_null_threshold)
      s$mu[l, ] = res$mu
      s$alpha[l, ] = res$alpha
      s$mu2[l, ] = res$mu2
      s$V[l,] = res$V #different from susieR
      s$lbf[l] = res$lbf_model
      s$KL[l] = -res$lbf_model + SER_posterior_e_loglik_rss(R,
                                                            Sigma, r, res$alpha * res$mu, res$alpha * res$mu2)
      s$Rz = s$Rz + R %*% (s$alpha[l, ] * s$mu[l, ])
    }
  s$Rz = unname(as.matrix(s$Rz))
  return(s)
}

assignInNamespace("update_each_effect_rss", update_each_effect, "susieR")

#' a susie fit object in order to initialize susie model. Revised from
#' susieR::init_finalize
init_finalize = function (s, X = NULL, Xr = NULL) {
  # different form susieR
  # if(length(s$V) == 1)
  #   s$V = rep(s$V, nrow(s$alpha))

  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")

  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")

  # check prior variance
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")

  # different from susieR
  # if (nrow(s$alpha) != length(s$V))
  #   stop("Input prior variance V must have length of nrow of alpha in ",
  #        "input object")

  # Update Xr.
  if (!missing(Xr))
    s$Xr = Xr
  if (!missing(X))
    s$Xr = susieR:::compute_Xb(X,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}

assignInNamespace("init_finalize", init_finalize ,"susieR")


#' a susie fit object in order to initialize susie rss model. Revised from
#' susieR::init_finalize_rss
init_finalize_rss <- function (s, R = NULL, Rz = NULL){
  # different from susieR
  # if (length(s$V) == 1)
  #   s$V = rep(s$V, nrow(s$alpha))
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("residual variance sigma2 must be positive (is your var(Y) zero?)")
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")

  # different from susieR
  # if (nrow(s$alpha) != length(s$V))
  #   stop("Input prior variance V must have length of nrow of alpha in ",
  #        "input object")

  if (!missing(Rz))
    s$Rz = Rz
  if (!missing(R))
    s$Rz = compute_Xb(R, colSums(s$mu * s$alpha))
  s$KL = rep(as.numeric(NA), nrow(s$alpha))
  s$lbf = rep(as.numeric(NA), nrow(s$alpha))
  class(s) = "susie"
  return(s)
}

assignInNamespace("init_finalize_rss", init_finalize ,"susieR")

#' susie get credible set. Revised from susieR::susie_get_cs
susie_get_cs <- function (res, X = NULL, Xcorr = NULL, coverage = 0.95, min_abs_corr = 0.5,
                          dedup = TRUE, squared = FALSE)
{
  if (!is.null(X) && !is.null(Xcorr)) {
    stop("Only one of X or Xcorr should be specified")
  }
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
    stop("Xcorr matrix must be symmetric")
  }
  if (inherits(res, "susie")) {
    null_index = res$null_index
    if (is.numeric(res$V))
      include_idx = rep(TRUE, nrow(res$alpha)) #different from susieR
    else include_idx = rep(TRUE, nrow(res$alpha))
  }
  else null_index = 0
  status = susieR:::in_CS(res$alpha, coverage)
  cs = lapply(1:nrow(status), function(i) which(status[i, ] !=
                                                  0))
  include_idx = include_idx * (lapply(cs, length) > 0)
  if (dedup)
    include_idx = include_idx * (!duplicated(cs))
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL, coverage = coverage))
  cs = cs[include_idx]
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) = paste0("L", which(include_idx))
    return(list(cs = cs, coverage = coverage))
  }
  else {
    purity = data.frame(do.call(rbind, lapply(1:length(cs),
                                              function(i) {
                                                if (null_index > 0 && null_index %in% cs[[i]])
                                                  c(-9, -9, -9)
                                                else susieR:::get_purity(cs[[i]], X, Xcorr, squared)
                                              })))
    if (squared)
      colnames(purity) = c("min.sq.corr", "mean.sq.corr",
                           "median.sq.corr")
    else colnames(purity) = c("min.abs.corr", "mean.abs.corr",
                              "median.abs.corr")
    threshold = ifelse(squared, min_abs_corr^2, min_abs_corr)
    is_pure = which(purity[, 1] >= threshold)
    if (length(is_pure) > 0) {
      cs = cs[is_pure]
      purity = purity[is_pure, ]
      row_names = paste0("L", which(include_idx)[is_pure])
      names(cs) = row_names
      rownames(purity) = row_names
      ordering = order(purity[, 1], decreasing = T)
      return(list(cs = cs[ordering], purity = purity[ordering,
                                                     ], cs_index = which(include_idx)[is_pure[ordering]],
                  coverage = coverage))
    }
    else {
      return(list(cs = NULL, coverage = coverage))
    }
  }
}

assignInNamespace("susie_get_cs", susie_get_cs ,"susieR")

#' susie get PIP. Revised from susieR::susie_get_pip
susie_get_pip <- function (res, prune_by_cs = FALSE, prior_tol = 1e-09)
{
  if (inherits(res, "susie")) {
    if (res$null_index > 0)
      res$alpha = res$alpha[, -res$null_index, drop = FALSE]
    if (is.numeric(res$V))
      include_idx = 1:nrow(res$alpha) # different from susieR, will affect value of small PIPs
    else include_idx = 1:nrow(res$alpha)
    if (!is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx, res$sets$cs_index)
    if (is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)
    if (length(include_idx) > 0) {
      res = res$alpha[include_idx, , drop = FALSE]
    }
    else {
      res = matrix(0, 1, ncol(res$alpha))
    }
  }
  return(as.vector(1 - apply(1 - res, 2, prod)))
}
assignInNamespace("susie_get_pip", susie_get_pip ,"susieR")



#' susie_rss. Revised from susieR::susie_rss
susie_rss <- function (z, R, maf = NULL, maf_thresh = 0, z_ld_weight = 0,
          L = 10, prior_variance = 50, residual_variance = NULL, r_tol = 1e-08,
          prior_weights = NULL, null_weight = NULL, estimate_residual_variance = TRUE,
          estimate_prior_variance = TRUE, estimate_prior_method = c("optim",
                                                                    "EM", "simple"), check_null_threshold = 0, prior_tol = 1e-09,
          max_iter = 100, s_init = list(), intercept_value = 0, coverage = 0.95,
          min_abs_corr = 0.5, tol = 0.001, verbose = FALSE, track_fit = FALSE,
          check_R = TRUE, check_z = TRUE)
{
  if (nrow(R) != length(z)) {
    stop(paste0("The dimension of correlation matrix (",
                nrow(R), " by ", ncol(R), ") does not agree with expected (",
                length(z), " by ", length(z), ")"))
  }
  if (!susieR:::is_symmetric_matrix(R)) {
    stop("R is not a symmetric matrix.")
  }
  if (!(is.double(R) & is.matrix(R)) & !inherits(R, "CsparseMatrix"))
    stop("Input R must be a double-precision matrix, or a sparse matrix.")
  if (!is.null(maf)) {
    if (length(maf) != length(z)) {
      stop(paste0("The length of maf does not agree with expected ",
                  length(z)))
    }
    id = which(maf > maf_thresh)
    R = R[id, id]
    z = z[id]
  }
  if (any(is.infinite(z))) {
    stop("z contains infinite value.")
  }
  if (any(is.na(R))) {
    stop("R matrix contains missing values.")
  }
  if (any(is.na(z))) {
    warning("NA values in z-scores are replaced with 0.")
    z[is.na(z)] = 0
  }
  if (z_ld_weight > 0) {
    R = susieR:::muffled_cov2cor((1 - z_ld_weight) * R + z_ld_weight *
                          tcrossprod(z))
    R = (R + t(R))/2
    check_z = FALSE
  }
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R) * (1 - null_weight),
                            ncol(R)), null_weight)
    else prior_weights = c(prior_weights * (1 - null_weight),
                           null_weight)
    R = cbind(rbind(R, 0), 0)
    z = c(z, 0)
  }
  if (!is.null(residual_variance) && (residual_variance >
                                      1 | residual_variance < 0)) {
    stop("Residual variance should be a scaler between 0 and 1.")
  }
  if (is.null(residual_variance)) {
    residual_variance = 1
  }
  p = ncol(R)
  attr(R, "eigen") = eigen(R, symmetric = TRUE)
  if (check_R && any(attr(R, "eigen")$values < -r_tol)) {
    stop(paste0("The correlation matrix (", nrow(R), " by ",
                ncol(R), "is not a positive semidefinite matrix. The smallest eigen value is ",
                min((attr(R, "eigen")$values), ". You can bypass this by \"check_R = FALSE\" which instead sets negative eigenvalues to 0 to allow for continued computations.")))
  }
  if (check_z) {
    proj = susieR:::check_projection(R, z)
    if (!proj$status) {
      warning("Input z does not lie in the space of non-zero eigenvectors of R. The result is thus not reliable.\n              Please refer to https://github.com/stephenslab/susieR/issues/91 for a possible solution.")
    }
    else {
      write("Input z is in space spanned by the non-zero eigenvectors of R.\n            You can safely set \"check_z = FALSE\" when you rerun the analysis, to save computation.",
            stderr())
    }
  }
  R = susieR:::set_R_attributes(R, r_tol)
  X = t(attr(R, "eigen")$vectors[, attr(R, "eigen")$values !=
                                   0]) * attr(R, "eigen")$values[attr(R, "eigen")$values !=
                                                                   0]^(0.5)
  Y = (t(attr(R, "eigen")$vectors[, attr(R, "eigen")$values !=
                                    0]) * attr(R, "eigen")$values[attr(R, "eigen")$values !=
                                                                    0]^(-0.5)) %*% z

  # different from susieR (as.numeric(var(Y)))
  s = susie(X, Y, L = L, scaled_prior_variance = prior_variance/as.numeric(var(Y)),
            residual_variance = residual_variance, prior_weights = prior_weights,
            null_weight = NULL, standardize = FALSE, intercept = FALSE,
            estimate_residual_variance = estimate_residual_variance,
            estimate_prior_variance = estimate_prior_variance, estimate_prior_method = estimate_prior_method,
            check_null_threshold = check_null_threshold, prior_tol = prior_tol,
            residual_variance_upperbound = 1, s_init = s_init, coverage = coverage,
            min_abs_corr = min_abs_corr, compute_univariate_zscore = FALSE,
            na.rm = FALSE, max_iter = max_iter, tol = tol, verbose = verbose,
            track_fit = track_fit)

  s$Rz = crossprod(X, s$Xr)
  s$fitted = s$Rz

  # different from susieR, this is a bug of susie_rss, make it consistent to susie
  if (!is.null(null_weight)) {
    s$null_index <- length(z)
    s$pip <- s$pip[1:(length(z)-1)]
  }

  return(s)
}

assignInNamespace("susie_rss", susie_rss ,"susieR")
