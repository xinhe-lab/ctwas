# Runs cTWAS fine-mapping for regions

Runs cTWAS fine-mapping for regions

## Usage

``` r
finemap_regions(
  region_data,
  LD_map,
  weights,
  L = 5,
  group_prior = NULL,
  group_prior_var = NULL,
  min_var = 2,
  min_gene = 1,
  null_method = c("ctwas", "susie", "none"),
  coverage = 0.95,
  min_abs_corr = 0.1,
  include_cs = TRUE,
  include_prior = FALSE,
  include_mu2 = FALSE,
  include_susie_alpha = TRUE,
  include_susie_result = FALSE,
  snps_only = FALSE,
  force_compute_cor = FALSE,
  save_cor = FALSE,
  cor_dir = NULL,
  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
  LD_loader_fun = NULL,
  snpinfo_loader_fun = NULL,
  ncore = 1,
  verbose = FALSE,
  logfile = NULL,
  ...
)
```

## Arguments

  - region\_data:
    
    region\_data to be finemapped

  - LD\_map:
    
    a data frame with filenames of LD matrices for the regions.

  - weights:
    
    a list of preprocessed weights.

  - L:
    
    the number of effects or a vector of number of effects for each
    region.

  - group\_prior:
    
    a vector of prior inclusion probabilities for different groups. If
    NULL, it will use uniform prior inclusion probabilities.

  - group\_prior\_var:
    
    a vector of prior variances for different groups. If NULL, it will
    set prior variance = 50 as the default in `susie_rss`.

  - min\_var:
    
    minimum number of variables (SNPs and genes) in a region.

  - min\_gene:
    
    minimum number of genes in a region.

  - null\_method:
    
    Method to compute null model, options: "ctwas", "susie" or "none".

  - coverage:
    
    A number between 0 and 1 specifying the “coverage” of the estimated
    confidence sets

  - min\_abs\_corr:
    
    Minimum absolute correlation allowed in a credible set.

  - include\_cs:
    
    If TRUE, add credible sets (CS) to fine-mapping results.

  - include\_prior:
    
    If TRUE, include priors in fine-mapping results.

  - include\_mu2:
    
    If TRUE, include estimated effect size variance (mu2) in
    fine-mapping results.

  - include\_susie\_alpha:
    
    If TRUE, include susie alpha matrix from fine-mapping results.

  - include\_susie\_result:
    
    If TRUE, include the "susie" result object in fine-mapping results.

  - snps\_only:
    
    If TRUE, use only SNPs in the region data.

  - force\_compute\_cor:
    
    If TRUE, force computing correlation (R) matrices

  - save\_cor:
    
    If TRUE, save correlation (R) matrices to `cor_dir`

  - cor\_dir:
    
    a string, the directory to store correlation (R) matrices

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
    
    The number of cores used to parallelize computation over regions

  - verbose:
    
    If TRUE, print detail messages

  - logfile:
    
    the log file, if NULL will print log info on screen

  - ...:
    
    Additional arguments of `susie_rss`.

## Value

a list with cTWAS fine-mapping results.
