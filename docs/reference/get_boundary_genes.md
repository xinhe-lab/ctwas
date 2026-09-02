# Get cross-boundary genes.

Get cross-boundary genes.

## Usage

``` r
get_boundary_genes(
  region_info,
  weights,
  gene_ids = NULL,
  mapping_table = NULL,
  show_mapping = TRUE,
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

  - mapping\_table:
    
    a data frame of mapping between molecular traits and genes, with
    required columns: "molecular\_id", "gene\_name".

  - show\_mapping:
    
    If TRUE, include the mapping between molecular traits and genes in
    the result.

  - ncore:
    
    The number of cores used to parallelize over genes.

## Value

a data frame of cross-boundary genes.
