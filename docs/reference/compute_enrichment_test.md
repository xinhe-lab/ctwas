# Computes enrichment (log-scale), standard error and p-value

Computes enrichment (log-scale), standard error and p-value

## Usage

``` r
compute_enrichment_test(
  group_prior,
  group_size,
  enrichment_test = c("fisher", "G", "z"),
  alternative = c("greater", "two.sided", "less"),
  include_test_result = FALSE
)
```

## Arguments

  - group\_prior:
    
    a vector of prior inclusion probabilities for different groups.

  - group\_size:
    
    a vector of number of variables in different groups.

  - enrichment\_test:
    
    Method to test enrichment, "G": G-test, "fisher": Fisher's exact
    test, "z": p-value computed from z-score of enrichment.

  - alternative:
    
    indicates the alternative hypothesis and must be one of "greater",
    "two.sided", or "less". Only used when `enrichment_test` is "fisher"
    or "z".

  - include\_test\_result:
    
    If TRUE, return the original test result.

## Value

Estimated enrichment, S.E. and p-value from G-test or Fisher's exact
test.
