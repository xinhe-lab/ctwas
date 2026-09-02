# Summarizes estimated parameters

Summarizes estimated parameters

## Usage

``` r
summarize_param(
  param,
  gwas_n,
  enrichment_test = c("fisher", "G", "z"),
  alternative = c("greater", "two.sided", "less")
)
```

## Arguments

  - param:
    
    a list of parameter estimation result from `est_param`

  - gwas\_n:
    
    the sample size of the GWAS summary statistics.

  - enrichment\_test:
    
    Method to test enrichment, "G": G-test, "fisher": Fisher's exact
    test, "z": p-value computed from z-score of enrichment.

  - alternative:
    
    indicates the alternative hypothesis and must be one of "greater",
    "two.sided", or "less". Only used when `enrichment_test` is "fisher"
    or "z".

## Value

a list of summarized parameters
