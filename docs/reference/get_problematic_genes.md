# Gets problematic genes from problematic SNPs

Gets problematic genes from problematic SNPs

## Usage

``` r
get_problematic_genes(
  problematic_snps,
  weights,
  finemap_res,
  pip_thresh = 0.5,
  z_thresh = NULL
)
```

## Arguments

  - problematic\_snps:
    
    a character vector of problematic SNP rsIDs

  - weights:
    
    a list of weights

  - finemap\_res:
    
    a data frame of cTWAS finemapping results

  - pip\_thresh:
    
    cutoff of gene PIP to select genes

  - z\_thresh:
    
    cutoff of abs(z-scores) to select genes

## Value

a vector of problematic genes
