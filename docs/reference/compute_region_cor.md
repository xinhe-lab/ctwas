# Computes correlation matrices for a single region

Computes correlation matrices for a single region

## Usage

``` r
compute_region_cor(sids, gids, R_snp, LD_sids, weights)
```

## Arguments

  - sids:
    
    SNP IDs in the region

  - gids:
    
    gene IDs in the region

  - R\_snp:
    
    LD (R) matrix

  - LD\_sids:
    
    SNP IDs for the rows and columns of the LD matrix

  - weights:
    
    a list of weights for all the genes

## Value

a list of correlation matrices (R\_snp, R\_snp\_gene and R\_gene)
