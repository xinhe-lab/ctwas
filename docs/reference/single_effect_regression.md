# Bayesian single-effect linear regression

These methods fit the regression model \\(y = Xb + e\\), where elements
of e are *i.i.d.* \\(N(0,s^2)\\), and b is a p-vector of effects to be
estimated. The assumption is that b has exactly one non-zero element,
with all elements equally likely to be non-zero. The prior on the
coefficient of the non-zero element is \\(N(0,V)\\).

## Usage

``` r
single_effect_regression(
  Y,
  X,
  V,
  residual_variance = 1,
  prior_weights = NULL,
  optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
  check_null_threshold = 0
)

single_effect_regression_rss(
  z,
  Sigma,
  V = 1,
  prior_weights = NULL,
  optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
  check_null_threshold = 0
)

single_effect_regression_ss(
  Xty,
  dXtX,
  V = 1,
  residual_variance = 1,
  prior_weights = NULL,
  optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
  check_null_threshold = 0
)
```

## Arguments

  - Y:
    
    An n-vector.

  - X:
    
    An n by p matrix of covariates.

  - V:
    
    A scalar giving the (initial) prior variance

  - residual\_variance:
    
    The residual variance.

  - prior\_weights:
    
    A p-vector of prior weights.

  - optimize\_V:
    
    The optimization method to use for fitting the prior variance.

  - check\_null\_threshold:
    
    Scalar specifying threshold on the log-scale to compare likelihood
    between current estimate and zero the null.

  - z:
    
    A p-vector of z scores.

  - Sigma:
    
    `residual_var*R + lambda*I`

  - Xty:
    
    A p-vector.

  - dXtX:
    
    A p-vector containing the diagonal elements of `crossprod(X)`.

## Value

A list with the following elements:

  - alpha:
    
    Vector of posterior inclusion probabilities; `alpha[i]` is posterior
    probability that the ith coefficient is non-zero.

  - mu:
    
    Vector of posterior means (conditional on inclusion).

  - mu2:
    
    Vector of posterior second moments (conditional on inclusion).

  - lbf:
    
    Vector of log-Bayes factors for each variable.

  - lbf\_model:
    
    Log-Bayes factor for the single effect regression.

`single_effect_regression` and `single_effect_regression_ss`
additionally output:

  - V:
    
    Prior variance (after optimization if `optimize_V != "none"`).

  - loglik:
    
    The log-likelihood, \\(\\log p(y | X, V)\\).

## Details

`single_effect_regression_ss` performs single-effect linear regression
with summary data, in which only the statistcs \\(X^Ty\\) and diagonal
elements of \\(X^TX\\) are provided to the method.

`single_effect_regression_rss` performs single-effect linear regression
with z scores. That is, this function fits the regression model \\(z =
R\*b + e\\), where e is \\(N(0,Sigma)\\), \\(Sigma = residual\_var\*R +
lambda\*I\\), and the b is a p-vector of effects to be estimated. The
assumption is that b has exactly one non-zero element, with all elements
equally likely to be non-zero. The prior on the non-zero element is
\\(N(0,V)\\). The required summary data are the p-vector `z` and the p
by p matrix `Sigma`. The summary statistics should come from the same
individuals.
