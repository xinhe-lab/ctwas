# Map finemapping result of molecular traits to genes.

Map finemapping result of molecular traits to genes.

## Usage

``` r
anno_finemap_res(
  finemap_res,
  snp_map,
  mapping_table,
  map_by = "molecular_id",
  add_gene_annot = TRUE,
  add_position = TRUE,
  use_gene_pos = c("mid", "start", "end"),
  drop_unmapped = TRUE
)
```

## Arguments

  - finemap\_res:
    
    a data frame of cTWAS finemapping results.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - mapping\_table:
    
    a data frame of mapping between molecular traits and genes, with
    required columns: "molecular\_id", "gene\_name".

  - map\_by:
    
    column name to be mapped by (default: "molecular\_id").

  - add\_gene\_annot:
    
    If TRUE, add annotations

  - add\_position:
    
    If TRUE, add positions

  - use\_gene\_pos:
    
    Use mid (midpoint), start, or end positions as gene positions.

  - drop\_unmapped:
    
    If TRUE, remove unmapped genes.

## Value

a data frame of annotated finemapping result.
