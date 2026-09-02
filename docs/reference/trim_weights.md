# Trim weights and filter SNPs by LD reference and GWAS SNPs.

Trim weights and filter SNPs by LD reference and GWAS SNPs.

## Usage

``` r
trim_weights(
  weights,
  snp_map,
  gwas_snp_ids,
  top_n_snps = NULL,
  min_abs_weight = 0,
  ncore = 1
)
```

## Arguments

  - weights:
    
    a list of pre-processed prediction weights.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - gwas\_snp\_ids:
    
    a vector of SNP IDs in GWAS summary statistics (z\_snp$id).

  - top\_n\_snps:
    
    an integer, only keeping the top n SNPs included in weight models.
    By default, keep all SNPs in weights.

  - min\_abs\_weight:
    
    a numeric value, only keeping the SNPs with abs(weight) \>
    `min_abs_weight`. By default, keep SNPs with abs(weight) \> 0.

  - ncore:
    
    The number of cores used to parallelize computation.

## Value

a list of trimmed weights.
