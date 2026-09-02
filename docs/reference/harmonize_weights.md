# Harmonize weight to match LD reference genotypes. Flip signs when reverse complement matches, remove ambiguous variants from the prediction models

Harmonize weight to match LD reference genotypes. Flip signs when
reverse complement matches, remove ambiguous variants from the
prediction models

## Usage

``` r
harmonize_weights(wgt.matrix, wgt.snpinfo, snp_info, drop_strand_ambig = TRUE)
```

## Arguments

  - wgt.matrix:
    
    a matrix of the weights

  - wgt.snpinfo:
    
    a data frame of the weight variants with columns "chrom", "id",
    "cm", "pos", "alt", "ref". "alt" is the effect allele.

  - snp\_info:
    
    a data frame of SNP info for the reference, with columns "chrom",
    "id", "pos", "alt", "ref".

  - drop\_strand\_ambig:
    
    If TRUE, remove strand ambiguous variants (A/T, G/C).

## Value

wgt.matrix and wgt.snpinfo with alleles flipped to match
