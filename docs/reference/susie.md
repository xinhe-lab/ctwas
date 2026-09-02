# Sum of Single Effects (SuSiE) Regression

Performs Bayesian multiple linear regression of Y on X; that is, this
function fits the regression model \\(Y = \\sum\_l X b\_{l=1}^L + e\\),
where elements of e are *i.i.d.* normal with zero mean and variance
`residual_variance`, and \\(\\sum\_{l=1}^L b\_l\\) is a vector of length
p representing the effects to be estimated. The “susie assumption” is
that each \\(b\_l\\) has exactly one non-zero element. The prior on the
non-zero element is normal with zero mean and variance `var(Y) *
scaled_prior_variance`. The model is fitted using the “Iterative
Bayesian Stepwise Selection” (IBSS) algorithm.

## Usage

``` r
susie(
  X,
  Y,
  L = min(10, ncol(X)),
  scaled_prior_variance = 0.2,
  residual_variance = NULL,
  prior_weights = NULL,
  null_weight = NULL,
  standardize = TRUE,
  intercept = TRUE,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = c("optim", "EM", "simple"),
  check_null_threshold = 0,
  prior_tol = 1e-09,
  residual_variance_upperbound = Inf,
  s_init = NULL,
  coverage = 0.95,
  min_abs_corr = 0.5,
  compute_univariate_zscore = FALSE,
  na.rm = FALSE,
  max_iter = 100,
  tol = 0.001,
  verbose = FALSE,
  track_fit = FALSE,
  residual_variance_lowerbound = var(drop(Y))/10000,
  refine = FALSE,
  warn_converge_fail = TRUE
)

susie_suff_stat(
  bhat,
  shat,
  R,
  n,
  var_y,
  XtX,
  Xty,
  yty,
  maf = NULL,
  maf_thresh = 0,
  L = 10,
  scaled_prior_variance = 0.2,
  residual_variance = NULL,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = c("optim", "EM", "simple"),
  check_null_threshold = 0,
  prior_tol = 1e-09,
  r_tol = 1e-08,
  prior_weights = NULL,
  null_weight = NULL,
  standardize = TRUE,
  max_iter = 100,
  s_init = NULL,
  intercept_value = 0,
  coverage = 0.95,
  min_abs_corr = 0.5,
  tol = 0.001,
  verbose = FALSE,
  track_fit = FALSE,
  check_input = FALSE,
  refine = FALSE,
  warn_converge_fail = TRUE
)
```

## Arguments

  - X:
    
    An n by p matrix of covariates.

  - Y:
    
    The observed responses, a vector of length n.

  - L:
    
    Number of components (nonzero coefficients) in the susie regression
    model. If L is larger than the number of covariates, p, L is set to
    p.

  - scaled\_prior\_variance:
    
    The scaled prior variance. This is either a scalar or a vector of
    length `L`. The prior variance of each non-zero element of b is set
    to `var(Y) * scaled_prior_variance`. If `estimate_prior_variance =
    TRUE`, this provides initial estimates of the prior variances.

  - residual\_variance:
    
    Variance of the residual. If `estimate_residual_variance = TRUE`,
    this value provides the initial estimate of the residual variance.
    By default, it is `var(Y)`.

  - prior\_weights:
    
    A vector of length p, in which each entry gives the prior
    probability that corresponding column of X has a nonzero effect on
    the outcome, Y.

  - null\_weight:
    
    Prior probability of no effect (a number between 0 and 1, and cannot
    be exactly 1).

  - standardize:
    
    If `standardize = TRUE`, standardize the columns of X (or XtX and
    Xty) to unit variance prior to fitting. Note that
    `scaled_prior_variance` specifies the prior on the coefficients of X
    *after* standardization (if it is performed). If you do not
    standardize, you may need to think more carefully about specifying
    `scaled_prior_variance`. Whatever your choice, the coefficients
    returned by `coef` are given for `X` on the original input scale.
    Any column of `X` that has zero variance is not standardized.

  - intercept:
    
    If `intercept = TRUE`, the intercept is fitted; it `intercept =
    FALSE`, the intercept is set to zero. Setting `intercept = FALSE` is
    generally not recommended.

  - estimate\_residual\_variance:
    
    If `estimate_residual_variance = TRUE`, the residual variance is
    estimated, using `residual_variance` as an initial value. If
    `estimate_residual_variance = FALSE`, the residual variance is fixed
    to the value supplied by `residual_variance`.

  - estimate\_prior\_variance:
    
    If `estimate_prior_variance = TRUE`, the prior variance is estimated
    (this is a separate parameter for each of the L effects). If
    provided, `scaled_prior_variance` is then used as an initial value
    for the optimization. When `estimate_prior_variance = FALSE`, the
    prior variance for each of the L effects is determined by the value
    supplied to `scaled_prior_variance`.

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

  - residual\_variance\_upperbound:
    
    Upper limit on the estimated residual variance. It is only relevant
    when `estimate_residual_variance = TRUE`.

  - s\_init:
    
    A previous susie fit with which to initialize.

  - coverage:
    
    A number between 0 and 1 specifying the “coverage” of the estimated
    confidence sets.

  - min\_abs\_corr:
    
    Minimum absolute correlation allowed in a credible set. The default,
    0.5, corresponds to a squared correlation of 0.25, which is a
    commonly used threshold for genotype data in genetic studies.

  - compute\_univariate\_zscore:
    
    If `compute_univariate_zscore = TRUE`, the univariate regression
    z-scores are outputted for each variable.

  - na.rm:
    
    Drop any missing values in Y from both X and Y.

  - max\_iter:
    
    Maximum number of IBSS iterations to perform.

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

  - residual\_variance\_lowerbound:
    
    Lower limit on the estimated residual variance. It is only relevant
    when `estimate_residual_variance = TRUE`.

  - refine:
    
    If `refine = TRUE`, we use a procedure to help SuSiE get out of
    local optimum.

  - warn\_converge\_fail:
    
    If `warn_converge_fail = TRUE`, prints a warning message when IBSS
    algorithm does not converge.

  - bhat:
    
    A p-vector of estimated effects.

  - shat:
    
    A p-vector of standard errors.

  - R:
    
    A p by p symmetric, positive semidefinite matrix. It can be
    \\(X'X\\), the covariance matrix \\(X'X/(n-1)\\), or a correlation
    matrix. It should be estimated from the same samples used to compute
    `bhat` and `shat`. Using an out-of-sample matrix may produce
    unreliable results.

  - n:
    
    The sample size.

  - var\_y:
    
    The sample variance of y, defined as \\(y'y/(n-1)\\). When the
    sample variance cannot be provided, the coefficients (returned from
    `coef`) are computed on the "standardized" X, y scale.

  - XtX:
    
    A p by p matrix \\(X'X\\) in which the columns of X are centered to
    have mean zero.

  - Xty:
    
    A p-vector \\(X'y\\) in which y and the columns of X are centered to
    have mean zero.

  - yty:
    
    A scalar \\(y'y\\) in which y is centered to have mean zero.

  - maf:
    
    Minor allele frequency; to be used along with `maf_thresh` to filter
    input summary statistics.

  - maf\_thresh:
    
    Variants having a minor allele frequency smaller than this threshold
    are not used.

  - r\_tol:
    
    Tolerance level for eigenvalue check of positive semidefinite matrix
    of R.

  - intercept\_value:
    
    The intercept. (The intercept cannot be estimated from centered
    summary data.) This setting will be used by `coef` to assign an
    intercept value, mainly for consistency with `susie`. Set to `NULL`
    if you want `coef` not to include an intercept term (and so only a
    p-vector is returned).

  - check\_input:
    
    If `check_input = TRUE`, `susie_suff_stat` performs additional
    checks on `XtX` and `Xty`. The checks are: (1) check that `XtX` is
    positive semidefinite; (2) check that `Xty` is in the space spanned
    by the non-zero eigenvectors of `XtX`.

## Value

A `"susie"` object with some or all of the following elements:

  - alpha:
    
    An L by p matrix of posterior inclusion probabilites.

  - mu:
    
    An L by p matrix of posterior means, conditional on inclusion.

  - mu2:
    
    An L by p matrix of posterior second moments, conditional on
    inclusion.

  - Xr:
    
    A vector of length n, equal to `X %*% colSums(alpha * mu)`.

  - lbf:
    
    log-Bayes Factor for each single effect.

  - lbf\_variable:
    
    log-Bayes Factor for each variable and single effect.

  - intercept:
    
    Intercept (fixed or estimated).

  - sigma2:
    
    Residual variance (fixed or estimated).

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

  - z:
    
    A vector of univariate z-scores.

  - niter:
    
    Number of IBSS iterations that were performed.

  - converged:
    
    `TRUE` or `FALSE` indicating whether the IBSS converged to a
    solution within the chosen tolerance level.

`susie_suff_stat` returns also outputs:

  - XtXr:
    
    A p-vector of `t(X)` times the fitted values, `X %*%
    colSums(alpha*mu)`.

## Details

`susie_suff_stat` performs sum of single-effect linear regression with
summary statistics. The required summary data are either: `bhat`,
`shat`, the p by p symmetric, positive semidefinite correlation (or
covariance) matrix `R`, the sample size `n`, and the variance of y; or
the p by p matrix \\(X'X\\), the p-vector \\(X'y\\), the sum of squares
\\(y'y\\), and the sample size `n`. The summary statistics should come
from the same individuals. Both the columns of X and the vector y should
be centered to have mean zero before computing these summary statistics;
you may also want to scale each column of X and y to have variance 1
(see examples).

## References

G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple new
approach to variable selection in regression, with application to
genetic fine-mapping. *Journal of the Royal Statistical Society, Series
B* [doi:10.1101/501114](https://doi.org/10.1101/501114) .

## See also

`susie_rss`
