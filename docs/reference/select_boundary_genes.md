# Gets boundary genes and selects high PIP boundary genes

Gets boundary genes and selects high PIP boundary genes

## Usage

``` r
select_boundary_genes(
  region_info,
  weights,
  gene_ids,
  finemap_res,
  susie_alpha_res,
  combine_PIPs = TRUE,
  mapping_table = NULL,
  pip_thresh = 0.5,
  filter_cs = FALSE,
  ncore = 1
)
```

## Arguments

  - region\_info:
    
    a data frame of region definitions.

  - weights:
    
    a list of preprocessed weights.

  - gene\_ids:
    
    a vector of selected gene IDs (z\_gene$id). If specified, limits to
    these genes. Default: use all genes in weights.

  - finemap\_res:
    
    a data frame of original finemapping result.

  - susie\_alpha\_res:
    
    a data frame of original susie alpha result.

  - combine\_PIPs:
    
    if TRUE, select boundary genes after combining gene PIPs.

  - mapping\_table:
    
    a data frame of mapping between molecular traits and genes, with
    required columns: "molecular\_id", "gene\_name".

  - pip\_thresh:
    
    PIP cutoff for selecting boundary genes to merge regions.

  - filter\_cs:
    
    If TRUE, only select boundary genes in credible sets for region
    merge.

  - ncore:
    
    The number of cores used to parallelize computation over regions

## Value

a list with boundary genes and selected boundary genes
