# Estimate L for all regions by running finemapping with uniform prior

Estimate L for all regions by running finemapping with uniform prior

## Usage

``` r
estimate_region_L(
  region_data,
  LD_map,
  weights,
  init_L = 5,
  min_abs_corr = 0.1,
  null_method = c("ctwas", "susie", "none"),
  null_weight = NULL,
  snps_only = FALSE,
  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
  LD_loader_fun = NULL,
  snpinfo_loader_fun = NULL,
  ncore = 1,
  verbose = FALSE,
  ...
)
```

## Arguments

  - region\_data:
    
    a list object indexing regions, variants and genes.

  - LD\_map:
    
    a data frame with filenames of LD matrices for each of the regions.

  - weights:
    
    a list of preprocessed weights.

  - init\_L:
    
    upper bound of the number of causal signals

  - min\_abs\_corr:
    
    Minimum absolute correlation allowed in a credible set.

  - null\_method:
    
    Method to compute null model, options: "ctwas", "susie" or "none".

  - null\_weight:
    
    Prior probability of no effect (a number between 0 and 1, and cannot
    be exactly 1). Only used when `null_method = "susie"`.

  - snps\_only:
    
    If TRUE, use only SNPs in the region data.

  - LD\_format:
    
    file format for LD matrix. If "custom", use a user defined
    `LD_loader_fun()` function to load LD matrix.

  - LD\_loader\_fun:
    
    a user defined function to load LD matrix when `LD_format =
    "custom"`.

  - snpinfo\_loader\_fun:
    
    a user defined function to load SNP information file, if SNP
    information files are not in standard cTWAS reference format.

  - ncore:
    
    The number of cores used to parallelize susie over regions

  - verbose:
    
    If TRUE, print detail messages

  - ...:
    
    Additional arguments of `susie_rss`.

## Value

a vector of estimated L for all regions
