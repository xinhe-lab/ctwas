# Loads weights in PredictDB format

Loads weights in PredictDB format

## Usage

``` r
load_predictdb_weights(
  weight_file,
  filter_protein_coding_genes = TRUE,
  load_predictdb_LD = TRUE
)
```

## Arguments

  - weight\_file:
    
    a string, pointing path to weights in PredictDB format.

  - filter\_protein\_coding\_genes:
    
    TRUE/FALSE. If TRUE, keep protein coding genes only.

  - load\_predictdb\_LD:
    
    TRUE/FALSE. If TRUE, load pre-computed LD among weight variants.
