# Diagnose LD mismatch using SuSiE RSS

Diagnose LD mismatch using SuSiE RSS

## Usage

``` r
diagnose_LD_mismatch_susie(
  region_ids,
  z_snp,
  LD_map,
  gwas_n,
  p_diff_thresh = 5e-06,
  plot = TRUE,
  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
  LD_loader_fun = NULL,
  snpinfo_loader_fun = NULL,
  ncore = 1,
  logfile = NULL
)
```

## Arguments

  - region\_ids:
    
    A vector of region IDs to run diagnosis

  - z\_snp:
    
    A data frame with two columns: "id", "A1", "A2", "z". giving the z
    scores for snps. "A1" is effect allele. "A2" is the other allele.

  - LD\_map:
    
    a data frame with filenames of LD matrices and SNP information for
    each of the regions.

  - gwas\_n:
    
    integer, GWAS sample size.

  - p\_diff\_thresh:
    
    numeric, p-value cutoff for identifying problematic SNPs with
    significant difference between observed z-scores and estimated
    values.

  - plot:
    
    If TRUE, plot observed z score vs the expected value.

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
    
    integer, number of cores for parallel computing.

  - logfile:
    
    the log file, if NULL will print log info on screen

## Value

a list of problematic SNPs, flipped SNPs, and test statistics from
susie's \`kriging\_rss\` function
