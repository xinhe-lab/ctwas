# @title Update each effect once.
# update each effect rss (revised from susieR)
# @param R a p by p symmetric and positive semidefinite correlation matrix.
# @param z a p vector
# @param s_init a list with elements sigma2, V, alpha, mu, Xr
# @param Sigma sigma2*R + lambda I
# @param estimate_prior_variance boolean indicating whether to estimate
#   prior variance
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
#
#' @importFrom Matrix diag
#'
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
