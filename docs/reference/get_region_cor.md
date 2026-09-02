# Gets correlation matrices for a single region.

Loads precomputed correlation matrices if available in `cor_dir`,
otherwise, computes correlation matrices if no precomputed correlation
matrices, or force\_compute\_cor = TRUE.

## Usage

``` r
get_region_cor(
  region_id,
  sids,
  gids,
  LD_map,
  weights,
  force_compute_cor = FALSE,
  save_cor = FALSE,
  cor_dir = NULL,
  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
  LD_loader_fun = NULL,
  snpinfo_loader_fun = NULL,
  verbose = FALSE
)
```

## Arguments

  - region\_id:
    
    a character string of region id.

  - sids:
    
    SNP IDs in the region

  - gids:
    
    gene IDs in the region

  - LD\_map:
    
    a data frame with filenames of LD matrices and SNP information for
    the regions.

  - weights:
    
    a list of preprocessed weights

  - force\_compute\_cor:
    
    If TRUE, force computing correlation (R) matrices

  - save\_cor:
    
    If TRUE, save correlation (R) matrices to `cor_dir`

  - cor\_dir:
    
    a string, the directory to store correlation (R) matrices

  - LD\_format:
    
    file format for LD matrix. If "custom", use a user defined
    `LD_loader_fun()` function to load LD matrix.

  - LD\_loader\_fun:
    
    a user defined function to load LD matrix when `LD_format =
    "custom"`.

  - snpinfo\_loader\_fun:
    
    a user defined function to load SNP information file, if SNP
    information files are not in standard cTWAS reference format.

  - verbose:
    
    If TRUE, print detail messages

## Value

correlation matrices (R\_snp, R\_snp\_gene and R\_gene)
