# Get correlations between CS, using variable with maximum PIP from each CS

Get correlations between CS, using variable with maximum PIP from each
CS

## Usage

``` r
get_cs_correlation(res, X = NULL, Xcorr = NULL, max = FALSE)
```

## Arguments

  - X:
    
    n by p matrix of values of the p variables (covariates) in n
    samples. When provided, correlation between variables will be
    computed and used to remove CSs whose minimum correlation among
    variables is smaller than `min_abs_corr`.

  - Xcorr:
    
    p by p matrix of correlations between variables (covariates). When
    provided, it will be used to remove CSs whose minimum correlation
    among variables is smaller than `min_abs_corr`.
