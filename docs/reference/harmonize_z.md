# Harmonizes z scores from GWAS to match LD reference genotypes. Flips signs when reverse complement matches, and removes strand ambiguous SNPs

Harmonizes z scores from GWAS to match LD reference genotypes. Flips
signs when reverse complement matches, and removes strand ambiguous SNPs

## Usage

``` r
harmonize_z(z_snp, snp_info, drop_strand_ambig = TRUE)
```

## Arguments

  - z\_snp:
    
    a data frame, with columns "id", "A1", "A2" and "z". Z scores for
    every SNP. "A1" is the effect allele.

  - snp\_info:
    
    a data frame of SNP info for the reference, with columns "chrom",
    "id", "pos", "alt", "ref".

  - drop\_strand\_ambig:
    
    If TRUE, remove strand ambiguous variants (A/T, G/C).

## Value

a data frame, z\_snp with the "z" columns flipped to match LD reference.
