# Updates cTWAS input data with merged region data

Updates cTWAS input data with merged region data

## Usage

``` r
update_merged_region_data(
  region_data,
  merged_region_data,
  region_info,
  merged_region_info,
  LD_map,
  merged_LD_map,
  snp_map,
  merged_snp_map,
  merged_region_id_map
)
```

## Arguments

  - region\_data:
    
    a list of original region data.

  - merged\_region\_data:
    
    a list of merged region data.

  - region\_info:
    
    a data frame of original region definitions.

  - merged\_region\_info:
    
    a data frame of original region definitions.

  - LD\_map:
    
    a data frame of original LD map.

  - merged\_LD\_map:
    
    a data frame of merged LD map.

  - snp\_map:
    
    a list of original SNP info.

  - merged\_snp\_map:
    
    a list of merged SNP info.

  - merged\_region\_id\_map:
    
    a data frame of new region IDs and original regions IDs.

## Value

a list with updated region\_data, region\_info, LD\_map, snp\_map.
