# Expands region\_data with all SNPs

Expands region\_data with all SNPs

## Usage

``` r
expand_region_data(region_data, snp_map, z_snp, maxSNP = Inf, ncore = 1)
```

## Arguments

  - region\_data:
    
    a list of assembled region data.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - maxSNP:
    
    Inf or integer. Maximum number of SNPs in a region. Default is Inf,
    no limit. This can be useful if there are many SNPs in a region and
    you don't have enough memory to run the program.

  - ncore:
    
    The number of cores used to parallelize susie over regions

## Value

updated region\_data with all SNPs
