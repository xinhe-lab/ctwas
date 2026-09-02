# Runs cTWAS post-processing procedure for region merging without LD

Runs cTWAS post-processing procedure for region merging without LD

## Usage

``` r
postprocess_region_merging_noLD(
  region_info,
  region_data,
  z_snp,
  z_gene,
  weights,
  snp_map,
  finemap_res,
  susie_alpha_res,
  combine_PIPs = TRUE,
  mapping_table = NULL,
  group_prior = NULL,
  group_prior_var = NULL,
  pip_thresh = 0.5,
  filter_cs = FALSE,
  maxSNP = Inf,
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
