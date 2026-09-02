# Inferences From Fitted SuSiE Model

These functions access basic properties or draw inferences from a fitted
susie model.

## Usage

``` r
susie_get_objective(res, last_only = TRUE, warning_tol = 1e-06)

susie_get_posterior_mean(res, prior_tol = 1e-09)

susie_get_posterior_sd(res, prior_tol = 1e-09)

susie_get_niter(res)

susie_get_prior_variance(res)

susie_get_residual_variance(res)

susie_get_lfsr(res)

susie_get_posterior_samples(susie_fit, num_samples)

susie_get_cs(
  res,
  X = NULL,
  Xcorr = NULL,
  coverage = 0.95,
  min_abs_corr = 0.5,
  dedup = TRUE,
  squared = FALSE
)

susie_get_pip(res, prune_by_cs = FALSE, prior_tol = 1e-09)
```

## Arguments

  - res:
    
    A susie fit, typically an output from `susie` or one of its
    variants. For `susie_get_pip` and `susie_get_cs`, this may instead
    be the posterior inclusion probability matrix, `alpha`.

  - last\_only:
    
    If `last_only = FALSE`, return the ELBO from all iterations;
    otherwise return the ELBO from the last iteration only.

  - warning\_tol:
    
    Warn if ELBO is decreasing by this tolerance level.

  - prior\_tol:
    
    Filter out effects having estimated prior variance smaller than this
    threshold.

  - susie\_fit:
    
    A susie fit, an output from `susie`.

  - num\_samples:
    
    The number of draws from the posterior distribution.

  - X:
    
    n by p matrix of values of the p variables (covariates) in n
    samples. When provided, correlation between variables will be
    computed and used to remove CSs whose minimum correlation among
    variables is smaller than `min_abs_corr`.

  - Xcorr:
    
    p by p matrix of correlations between variables (covariates). When
    provided, it will be used to remove CSs whose minimum correlation
    among variables is smaller than `min_abs_corr`.

  - coverage:
    
    A number between 0 and 1 specifying desired coverage of each CS.

  - min\_abs\_corr:
    
    A "purity" threshold for the CS. Any CS that contains a pair of
    variables with correlation less than this threshold will be filtered
    out and not reported.

  - dedup:
    
    If `dedup = TRUE`, remove duplicate CSs.

  - squared:
    
    If `squared = TRUE`, report min, mean and median of squared
    correlation instead of the absolute correlation.

  - prune\_by\_cs:
    
    Whether or not to ignore single effects not in a reported CS when
    calculating PIP.

## Value

`susie_get_objective` returns the evidence lower bound (ELBO) achieved
by the fitted susie model and, optionally, at each iteration of the IBSS
fitting procedure.

`susie_get_residual_variance` returns the (estimated or fixed) residual
variance parameter.

`susie_get_prior_variance` returns the (estimated or fixed) prior
variance parameters.

`susie_get_posterior_mean` returns the posterior mean for the regression
coefficients of the fitted susie model.

`susie_get_posterior_sd` returns the posterior standard deviation for
coefficients of the fitted susie model.

`susie_get_niter` returns the number of model fitting iterations
performed.

`susie_get_pip` returns a vector containing the posterior inclusion
probabilities (PIPs) for all variables.

`susie_get_lfsr` returns a vector containing the average lfsr across
variables for each single-effect, weighted by the posterior inclusion
probability (alpha).

`susie_get_posterior_samples` returns a list containing the effect sizes
samples and causal status with two components: `b`, an `num_variables` x
`num_samples` matrix of effect sizes; `gamma`, an `num_variables` x
`num_samples` matrix of causal status random draws.

`susie_get_cs` returns credible sets (CSs) from a susie fit, as well as
summaries of correlation among the variables included in each CS. If
desired, one can filter out CSs that do not meet a specified “purity”
threshold; to do this, either `X` or `Xcorr` must be supplied. It
returns a list with the following elements:

  - cs:
    
    A list in which each list element is a vector containing the indices
    of the variables in the CS.

  - coverage:
    
    The nominal coverage specified for each CS.

  - purity:
    
    If `X` or `Xcorr` iis provided), the purity of each CS.

  - cs\_index:
    
    If `X` or `Xcorr` is provided) the index (number between 1 and L) of
    each reported CS in the supplied susie fit.
