# Preprocess GWAS z-scores, harmonize GWAS z-scores with LD reference

Preprocess GWAS z-scores, harmonize GWAS z-scores with LD reference

## Usage

``` r
preprocess_z_snp(
  z_snp,
  snp_map,
  drop_multiallelic = TRUE,
  drop_strand_ambig = TRUE,
  varID_converter_fun = NULL,
  logfile = NULL
)
```

## Arguments

  - z\_snp:
    
    A data frame with two columns: "id", "A1", "A2", "z". giving the z
    scores for snps. "A1" is effect allele. "A2" is the other allele.

  - snp\_map:
    
    a list of SNP-to-region map for the reference.

  - drop\_multiallelic:
    
    If TRUE, multiallelic variants will be dropped from the summary
    statistics.

  - drop\_strand\_ambig:
    
    If TRUE remove strand ambiguous variants (A/T, G/C).

  - varID\_converter\_fun:
    
    a user defined function to convert GWAS variant IDs to the reference
    variant format.

  - logfile:
    
    The log filename. If NULL, print log info on screen.

## Value

a data frame preprocessed GWAS z-scores
