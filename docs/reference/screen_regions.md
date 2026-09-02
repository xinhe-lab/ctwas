# Screens regions with strong signals

Screens regions with strong signals

## Usage

``` r
screen_regions(
  region_data,
  group_prior = NULL,
  group_prior_var = NULL,
  min_var = 2,
  min_gene = 1,
  min_nonSNP_PIP = 0.5,
  min_snp_pval = 5e-08,
  null_method = c("ctwas", "susie", "none"),
  ncore = 1,
  logfile = NULL,
  verbose = FALSE
)
```

## Arguments

  - region\_data:
    
    a list object indexing regions, variants and genes.

  - group\_prior:
    
    a vector of prior inclusion probabilities for different groups.

  - group\_prior\_var:
    
    a vector of prior variances for different groups.

  - min\_var:
    
    minimum number of variables (SNPs and genes) in a region.

  - min\_gene:
    
    minimum number of genes in a region.

  - min\_nonSNP\_PIP:
    
    If screening by non-SNP PIPs, regions with total non-SNP PIP \>=
    `min_nonSNP_PIP` will be selected to run finemapping using full
    SNPs.

  - min\_snp\_pval:
    
    Select regions with minimum SNP p-values \< `min_snp_pval`.

  - null\_method:
    
    Method to compute null model, options: "ctwas", "susie" or "none".

  - ncore:
    
    The number of cores used to parallelize susie over regions.

  - logfile:
    
    The log filename. If NULL, print log info on screen.

  - verbose:
    
    If TRUE, print detail messages.

## Value

a list, containing a data frame of selected region data, and a data
frame of screening summary for all regions.
