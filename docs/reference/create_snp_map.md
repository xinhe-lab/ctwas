# Map SNPs to regions using all the variants in the LD reference.

Map SNPs to regions using all the variants in the LD reference.

## Usage

``` r
create_snp_map(region_info, ref_snp_info, ncore = 1)
```

## Arguments

  - region\_info:
    
    a data frame of region definitions, with columns: "chrom", "start",
    "stop", "region\_id".

  - ref\_snp\_info:
    
    a data frame of all variant info in the reference.

  - ncore:
    
    The number of cores used to parallelize over regions
