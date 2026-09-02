# Map molecular traits to genes in susie alpha result.

Map molecular traits to genes in susie alpha result.

## Usage

``` r
anno_susie_alpha_res(
  susie_alpha_res,
  mapping_table,
  map_by = "molecular_id",
  drop_unmapped = TRUE
)
```

## Arguments

  - susie\_alpha\_res:
    
    a data frame of susie alpha result.

  - mapping\_table:
    
    a data frame of mapping between molecular traits and genes, with
    required columns: "molecular\_id", "gene\_name".

  - map\_by:
    
    column name to be mapped by (default: "molecular\_id").

  - drop\_unmapped:
    
    If TRUE, remove unmapped genes.

## Value

a data frame of annotated finemapping result.
