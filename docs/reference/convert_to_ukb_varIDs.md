# Convert variant IDs from Open GWAS format or PredictDB weight format

Convert variant IDs from Open GWAS format or PredictDB weight format

## Usage

``` r
convert_to_ukb_varIDs(varIDs, ref_format = "%s:%s_%s_%s")
```

## Arguments

  - varIDs:
    
    a vector of variant IDs from Open GWAS format ("chr\_pos\_ref\_alt")
    or PredictDB weight format (“chr\_pos\_ref\_alt\_build” or
    “chr:pos\_ref\_alt\_build”)

  - ref\_format:
    
    variant ID format in the LD reference.

## Value

a vector of variant IDs in the format of the LD reference.
