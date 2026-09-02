# Add SNP and gene positions to cTWAS finemapping result.

Add SNP and gene positions to cTWAS finemapping result.

## Usage

``` r
add_pos(
  finemap_res,
  snp_map,
  mapping_table,
  map_by = "molecular_id",
  use_gene_pos = c("mid", "start", "end")
)
```

## Arguments

  - finemap\_res:
    
    a data frame of cTWAS finemapping result

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - mapping\_table:
    
    a data frame of gene annotations, with required columns:
    "gene\_name", "chrom", "pos" (or "start" and "end").

  - map\_by:
    
    column name to be mapped by (default: "gene\_name").

  - use\_gene\_pos:
    
    Use mid (midpoint), start, or end positions as gene positions.

## Value

a data frame with chromosomes and positions added to fine-mapping result
