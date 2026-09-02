# Compute Distribution of z-scores of Variant j Given Other z-scores, and Detect Possible Allele Switch Issue

Under the null, the rss model with regularized LD matrix is \\(z|R,s \~
N(0, (1-s)R + s I))\\). We use a mixture of normals to model the
conditional distribution of z\_j given other z scores, \\(z\_j |
z\_{-j}, R, s \~ \\sum\_{k=1}^{K} \\pi\_k N(-\\Omega\_{j,-j}
z\_{-j}/\\Omega\_{jj}, \\sigma\_{k}^2/\\Omega\_{jj})\\), \\(\\Omega =
((1-s)R + sI)^{-1}\\), \\(\\sigma\_1, ..., \\sigma\_k\\) is a grid of
fixed positive numbers. We estimate the mixture weights \\(\\pi\\) We
detect the possible allele switch issue using likelihood ratio for each
variant.

## Usage

``` r
kriging_rss(
  z,
  R,
  n,
  r_tol = 1e-08,
  s = estimate_s_rss(z, R, n, r_tol, method = "null-mle"),
  plot = FALSE
)
```

## Arguments

  - z:
    
    A p-vector of z scores.

  - R:
    
    A p by p symmetric, positive semidefinite correlation matrix.

  - n:
    
    The sample size. (Optional, but highly recommended.)

  - r\_tol:
    
    Tolerance level for eigenvalue check of positive semidefinite matrix
    of R.

  - s:
    
    an estimated s from `estimate_s_rss`

  - plot:
    
    If TRUE, plot observed z score vs the expected value.

## Value

a list containing a ggplot2 plot object and a table. The plot compares
observed z score vs the expected value. The possible allele switched
variants are labeled as red points (log LR \> 2 and abs(z) \> 2). The
table summarizes the conditional distribution for each variant and the
likelihood ratio test. The table has the following columns: the observed
z scores, the conditional expectation, the conditional variance, the
standardized differences between the observed z score and expected
value, the log likelihood ratio statistics.
