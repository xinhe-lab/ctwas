# Map regions for each gene.

Map regions for each gene.

## Usage

``` r
map_gene_regions(gene_info, region_info, ncore = 1)
```

## Arguments

  - gene\_info:
    
    a data frame of gene info.

  - region\_info:
    
    a data frame of region definitions.

  - ncore:
    
    The number of cores used to parallelize over genes.

## Value

a data frame of gene info with regions overlapping each gene
