# Computes z-scores of molecular traits

Computes z-scores of molecular traits

## Usage

``` r
compute_gene_z(z_snp, weights, ncore = 1, logfile = NULL)
```

## Arguments

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - weights:
    
    a list of preprocessed weights.

  - ncore:
    
    The number of cores used to parallelize computation over weights.

  - logfile:
    
    The log filename. If NULL, print log info on screen.

## Value

a data frame of z-scores of molecular traits
