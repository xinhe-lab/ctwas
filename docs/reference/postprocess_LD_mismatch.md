# Runs cTWAS post-processing procedure for merging regions

Runs cTWAS post-processing procedure for merging regions

## Usage

``` r
postprocess_LD_mismatch(
  region_data,
  z_snp,
  weights,
  LD_map,
  snp_map,
  gwas_n,
  finemap_res,
  susie_alpha_res,
  group_prior = NULL,
  group_prior_var = NULL,
  min_nonSNP_PIP = 0.8,
  filter_cs = FALSE,
  p_diff_thresh = 5e-06,
  pip_thresh = 0.5,
  z_thresh = NULL,
  plot = TRUE,
  maxSNP = Inf,
  ncore = 1,
  verbose = FALSE,
  logfile = NULL,
  ...
)
```

## Arguments

  - region\_data:
    
    region\_data to be fine-mapped.

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - weights:
    
    a list of preprocessed weights.

  - LD\_map:
    
    a data frame with filenames of LD matrices for the regions.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - gwas\_n:
    
    integer, GWAS sample size.

  - finemap\_res:
    
    a data frame of original finemapping result.

  - susie\_alpha\_res:
    
    a data frame of original susie alpha result.

  - group\_prior:
    
    a vector of prior inclusion probabilities for different groups. If
    NULL, it will use uniform prior inclusion probabilities.

  - group\_prior\_var:
    
    a vector of prior variances for different groups. If NULL, it will
    set prior variance = 50 as the default in `susie_rss`.

  - min\_nonSNP\_PIP:
    
    regions with total non-SNP PIP \>= `min_nonSNP_PIP` will be selected
    to run LD mismatch diagnosis.

  - filter\_cs:
    
    If TRUE, computes region non-SNP PIP only in credible sets.

  - p\_diff\_thresh:
    
    p-value cutoff for identifying problematic SNPs with significant
    difference between observed z-scores and estimated values.

  - pip\_thresh:
    
    cutoff of gene PIP to select problematic genes.

  - z\_thresh:
    
    cutoff of abs(z-scores) to select problematic genes.

  - plot:
    
    If TRUE, plot observed z-score vs the expected value.

  - maxSNP:
    
    Inf or integer. Maximum number of SNPs in a region. Default is Inf,
    no limit. This can be useful if there are many SNPs in a region and
    you don't have enough memory to run the program.

  - ncore:
    
    The number of cores used to parallelize computation over regions

  - verbose:
    
    If TRUE, print detail messages.

  - logfile:
    
    the log file, if NULL will print log info on screen

  - ...:
    
    Additional arguments of `finemap_regions`.

## Value

a list with region merge results.
