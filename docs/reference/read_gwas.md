# Read GWAS summary statistics

Read GWAS summary statistics

## Usage

``` r
read_gwas(
  gwas,
  id = "rsid",
  A1 = "ALT",
  A2 = "REF",
  z = "Z",
  beta = "ES",
  se = "SE"
)
```

## Arguments

  - gwas:
    
    A data frame of GWAS summary statistics.

  - id:
    
    Column name of the variant IDs.

  - A1:
    
    Column name of the alternative alleles.

  - A2:
    
    Column name of the reference alleles.

  - z:
    
    Column name of z-scores.

  - beta:
    
    Column name of effect sizes. If z is not available, compute z from
    beta and se.

  - se:
    
    Column name of the standard errors.

## Value

A data frame of processed GWAS summary statistics.
