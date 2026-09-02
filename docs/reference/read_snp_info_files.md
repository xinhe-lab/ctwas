# Read multiple single SNP info files as a data frame

Read multiple single SNP info files as a data frame

## Usage

``` r
read_snp_info_files(files, snpinfo_loader_fun = NULL)
```

## Arguments

  - files:
    
    a vector of file names for SNP info in all regions

  - snpinfo\_loader\_fun:
    
    a user defined function to load SNP information file, if SNP
    information files are not in standard cTWAS reference format.

## Value

a data frame of SNP info from all files
