# convert z-scores to p-values (two-sided test)

convert z-scores to p-values (two-sided test)

## Usage

``` r
z2p(z, neg_log10_p = FALSE)
```

## Arguments

  - z:
    
    a vector of z-scores

  - neg\_log10\_p:
    
    If TRUE, returns -log10(p-values). Otherwise, returns p-values.

## Value

a vector of p-values, or -log10(p-values) if neg\_log10\_p is TRUE.
