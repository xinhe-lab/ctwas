# Read a single SNP info file as a data frame

Read a single SNP info file as a data frame

## Usage

``` r
read_snp_info_file(file, snpinfo_loader_fun = NULL)
```

## Arguments

  - file:
    
    SNP info file name.

  - snpinfo\_loader\_fun:
    
    a user defined function to load SNP information file, if SNP
    information files are not in standard cTWAS reference format.

## Value

a data frame of SNP info
