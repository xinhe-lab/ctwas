# Sum of Single Effects (SuSiE) Regression using summary statistics

`susie_rss` performs sum of single-effect linear regression with z
scores; all posterior calculations are for z-scores. This function fits
the regression model \\(z = \\sum\_l R\*b\_l + e\\), where e is
\\(N(0,R)\\) and \\(\\sum\_l b\_l\\) is a p-vector of effects to be
estimated. The required summary data are the p by p correlation matrix,
`R`, and the p-vector `z`.

## Usage

``` r
susie_rss(
  z,
  R,
  maf = NULL,
  maf_thresh = 0,
  z_ld_weight = 0,
  L = 10,
  prior_variance = 50,
  residual_variance = NULL,
  prior_weights = NULL,
  null_weight = NULL,
  estimate_residual_variance = FALSE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = c("optim", "EM", "simple"),
  check_null_threshold = 0,
  prior_tol = 1e-09,
  max_iter = 100,
  s_init = NULL,
  intercept_value = 0,
  coverage = 0.95,
  min_abs_corr = 0.5,
  tol = 0.001,
  verbose = FALSE,
  track_fit = FALSE,
  check_R = FALSE,
  r_tol = 1e-08,
  refine = FALSE,
  warn_converge_fail = TRUE
)
```

## Arguments

  - z:
    
    A p-vector of z scores.

  - R:
    
    A p by p symmetric, positive semidefinite correlation matrix.

  - maf:
    
    Minor allele frequency; to be used along with `maf_thresh` to filter
    input summary statistics.

  - maf\_thresh:
    
    Variants having a minor allele frequency smaller than this threshold
    are not used.

  - z\_ld\_weight:
    
    This feature is not recommended. The weights assigned to the z
    scores in the LD matrix. If `z_ld_weight > 0`, the LD matrix used in
    the model is `cov2cor((1-w)*R + w*tcrossprod(z))`, where `w =
    z_ld_weight`.

  - L:
    
    Number of components (nonzero coefficients) in the susie regression
    model. If L is larger than the number of covariates, p, L is set to
    p.

  - prior\_variance:
    
    The prior variance. It is either a scalar or a vector of length L.

  - residual\_variance:
    
    Variance of the residual. If it is not specified, we set it to 1.

  - prior\_weights:
    
    A vector of length p, in which each entry gives the prior
    probability that SNP j has non-zero effect.

  - null\_weight:
    
    Prior probability of no effect (a number between 0 and 1, and cannot
    be exactly 1).

  - estimate\_residual\_variance:
    
    The residual variance is fixed to the value supplied by
    `residual_variance`. We don't estimate residual variance in
    susie\_rss.

  - estimate\_prior\_variance:
    
    If `estimate_prior_variance = TRUE`, the prior variance is estimated
    (this is a separate parameter for each of the L effects). If
    provided, `prior_variance` is then used as an initial value for the
    optimization. When `estimate_prior_variance = FALSE`, the prior
    variance for each of the L effects is determined by the value
    supplied to `prior_variance`.

  - estimate\_prior\_method:
    
    The method used for estimating prior variance. When
    `estimate_prior_method = "simple"` is used, the likelihood at the
    specified prior variance is compared to the likelihood at a variance
    of zero, and the setting with the larger likelihood is retained.

  - check\_null\_threshold:
    
    When the prior variance is estimated, compare the estimate with the
    null, and set the prior variance to zero unless the log-likelihood
    using the estimate is larger by this threshold amount. For example,
    if you set `check_null_threshold = 0.1`, this will "nudge" the
    estimate towards zero when the difference in log-likelihoods is
    small. A note of caution that setting this to a value greater than
    zero may lead the IBSS fitting procedure to occasionally decrease
    the ELBO.

  - prior\_tol:
    
    When the prior variance is estimated, compare the estimated value to
    `prior_tol` at the end of the computation, and exclude a single
    effect from PIP computation if the estimated prior variance is
    smaller than this tolerance value.

  - max\_iter:
    
    Maximum number of IBSS iterations to perform.

  - s\_init:
    
    A previous susie fit with which to initialize.

  - intercept\_value:
    
    The intercept. (The intercept cannot be estimated from centered
    summary data.) This setting will be used by `coef` to assign an
    intercept value, mainly for consistency with `susie`. Set to `NULL`
    if you want `coef` not to include an intercept term (and so only a
    p-vector is returned).

  - coverage:
    
    A number between 0 and 1 specifying the “coverage” of the estimated
    confidence sets.

  - min\_abs\_corr:
    
    Minimum absolute correlation allowed in a credible set. The default,
    0.5, corresponds to a squared correlation of 0.25, which is a
    commonly used threshold for genotype data in genetic studies.

  - tol:
    
    A small, non-negative number specifying the convergence tolerance
    for the IBSS fitting procedure. The fitting procedure will halt when
    the difference in the variational lower bound, or “ELBO” (the
    objective function to be maximized), is less than `tol`.

  - verbose:
    
    If `verbose = TRUE`, the algorithm's progress, and a summary of the
    optimization settings, are printed to the console.

  - track\_fit:
    
    If `track_fit = TRUE`, `trace` is also returned containing detailed
    information about the estimates at each iteration of the IBSS
    fitting procedure.

  - check\_R:
    
    If `check_R = TRUE`, check that `R` is positive semidefinite.

  - r\_tol:
    
    Tolerance level for eigenvalue check of positive semidefinite matrix
    of R.

  - refine:
    
    If `refine = TRUE`, we use a procedure to help SuSiE get out of
    local optimum.

  - warn\_converge\_fail:
    
    If `warn_converge_fail = TRUE`, prints a warning message when IBSS
    algorithm does not converge.

## Value

A `"susie"` object with some or all of the following elements:

  - alpha:
    
    An L by p matrix of posterior inclusion probabilites.

  - mu:
    
    An L by p matrix of posterior means, conditional on inclusion.

  - mu2:
    
    An L by p matrix of posterior second moments, conditional on
    inclusion.

  - lbf:
    
    log-Bayes Factor for each single effect.

  - lbf\_variable:
    
    log-Bayes Factor for each variable and single effect.

  - intercept:
    
    Fixed Intercept.

  - sigma2:
    
    Fixed Residual variance.

  - V:
    
    Prior variance of the non-zero elements of b, equal to
    `scaled_prior_variance * var(Y)`.

  - elbo:
    
    The value of the variational lower bound, or “ELBO” (objective
    function to be maximized), achieved at each iteration of the IBSS
    fitting procedure.

  - fitted:
    
    Vector of length n containing the fitted values of the outcome.

  - sets:
    
    Credible sets estimated from model fit; see `susie_get_cs` for
    details.

  - pip:
    
    A vector of length p giving the (marginal) posterior inclusion
    probabilities for all p covariates.

  - niter:
    
    Number of IBSS iterations that were performed.

  - converged:
    
    `TRUE` or `FALSE` indicating whether the IBSS converged to a
    solution within the chosen tolerance level.

  - Rr:
    
    An p-vector of `t(X)` times fitted values, `X %*%
    colSums(alpha*mu)`.
