# Map SNPs to regions using region meta table.

Map SNPs to regions using region meta table.

## Usage

``` r
create_snp_LD_map(region_metatable, snpinfo_loader_fun = NULL)
```

## Arguments

  - region\_metatable:
    
    a data frame of region meta table, with columns: "chrom", "start",
    "stop", "region\_id", "LD\_file", "SNP\_file".

  - snpinfo\_loader\_fun:
    
    IF not NULL, use custom loader function to read SNP files.
