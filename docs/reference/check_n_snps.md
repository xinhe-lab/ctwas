# Check the numbers of SNPs in snp\_map, z\_snp and weights

Check the numbers of SNPs in snp\_map, z\_snp and weights

## Usage

``` r
check_n_snps(snp_map, z_snp, weights)
```

## Arguments

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - weights:
    
    a list of pre-processed prediction weights.

## Value

a list of numbers of SNPs in snp\_map, z\_snp and weights
