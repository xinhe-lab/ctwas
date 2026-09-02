# Updates cTWAS finemapping result for merged regions

Updates cTWAS finemapping result for merged regions

## Usage

``` r
update_merged_region_finemap_res(
  finemap_res,
  susie_alpha_res,
  merged_region_finemap_res,
  merged_region_susie_alpha_res,
  merged_region_id_map
)
```

## Arguments

  - finemap\_res:
    
    a data frame of original finemapping result.

  - susie\_alpha\_res:
    
    a data frame of original susie alpha result.

  - merged\_region\_finemap\_res:
    
    a data frame of finemapping result for merged regions.

  - merged\_region\_susie\_alpha\_res:
    
    a data frame of susie alpha result for merged regions.

  - merged\_region\_id\_map:
    
    a data frame of new region IDs and original regions IDs.

## Value

a list with updated cTWAS finemapping result.
