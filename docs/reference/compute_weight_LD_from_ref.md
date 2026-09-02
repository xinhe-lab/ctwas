# Computes LD for weight variants using reference LD

Computes LD for weight variants using reference LD

## Usage

``` r
compute_weight_LD_from_ref(
  weights,
  region_info,
  LD_map,
  snp_map = NULL,
  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
  LD_loader_fun = NULL,
  snpinfo_loader_fun = NULL,
  ncore = 1
)
```

## Arguments

  - weights:
    
    a list of preprocessed weights.

  - region\_info:
    
    a data frame of region definitions.

  - LD\_map:
    
    a data frame with filenames of LD matrices and SNP information for
    the regions. Required when `load_predictdb_LD = FALSE`.

  - snp\_map:
    
    a list of SNP-to-region map for the reference. If NUll, it will
    reads SNP info from the "SNP\_file" column of LD\_map.

  - LD\_format:
    
    file format for LD matrix. If "custom", use a user defined
    `LD_loader_fun()` function to load LD matrix.

  - LD\_loader\_fun:
    
    a user defined function to load LD matrix when `LD_format =
    "custom"`.

  - snpinfo\_loader\_fun:
    
    a user defined function to load SNP information file, if SNP
    information files are not in standard cTWAS reference format.

  - ncore:
    
    The number of cores used to parallelize computation.

## Value

a list of processed weights, with LD of weight variants included.
