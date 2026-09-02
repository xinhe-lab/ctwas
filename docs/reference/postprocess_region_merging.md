# Runs cTWAS post-processing procedure for region merging

Runs cTWAS post-processing procedure for region merging

## Usage

``` r
postprocess_region_merging(
  region_info,
  region_data,
  z_snp,
  z_gene,
  weights,
  LD_map,
  snp_map,
  finemap_res,
  susie_alpha_res,
  combine_PIPs = TRUE,
  mapping_table = NULL,
  L = 5,
  group_prior = NULL,
  group_prior_var = NULL,
  pip_thresh = 0.5,
  filter_cs = FALSE,
  maxSNP = Inf,
  save_cor = FALSE,
  cor_dir = NULL,
  ncore = 1,
  verbose = FALSE,
  logfile = NULL,
  ...
)
```

## Arguments

  - region\_info:
    
    a data frame of region definitions.

  - region\_data:
    
    region\_data to be fine-mapped.

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - z\_gene:
    
    A data frame with columns: "id", "z", giving the z-scores for genes.

  - weights:
    
    a list of preprocessed weights.

  - LD\_map:
    
    a data frame with filenames of LD matrices for the regions.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - finemap\_res:
    
    a data frame of original finemapping result.

  - susie\_alpha\_res:
    
    a data frame of original susie alpha result.

  - combine\_PIPs:
    
    if TRUE, select boundary genes after combining gene PIPs.

  - mapping\_table:
    
    a data frame of mapping between molecular traits and genes, with
    required columns: "molecular\_id", "gene\_name".

  - L:
    
    the number of effects or a vector of number of effects for each
    region.

  - group\_prior:
    
    a vector of prior inclusion probabilities for different groups. If
    NULL, it will use uniform prior inclusion probabilities.

  - group\_prior\_var:
    
    a vector of prior variances for different groups. If NULL, it will
    set prior variance = 50 as the default in `susie_rss`.

  - pip\_thresh:
    
    PIP cutoff for selecting boundary genes to merge regions.

  - filter\_cs:
    
    If TRUE, only select boundary genes in credible sets for region
    merge.

  - maxSNP:
    
    Inf or integer. Maximum number of SNPs in a region. Default is Inf,
    no limit. This can be useful if there are many SNPs in a region and
    you don't have enough memory to run the program.

  - save\_cor:
    
    If TRUE, save correlation (R) matrices to `cor_dir`

  - cor\_dir:
    
    a string, the directory to store correlation (R) matrices

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
