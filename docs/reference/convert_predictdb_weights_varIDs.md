# Convert rsIDs and varIDs in predictDB weights to the reference variant format.

Convert rsIDs and varIDs in predictDB weights to the reference variant
format.

## Usage

``` r
convert_predictdb_weights_varIDs(
  weight_file,
  load_predictdb_LD = TRUE,
  varID_converter_fun = NULL,
  rsID_converter_fun = NULL,
  outputdir = getwd(),
  outname
)
```

## Arguments

  - weight\_file:
    
    filename of the PredictDB weights '.db' file.

  - load\_predictdb\_LD:
    
    If TRUE, load pre-computed LD among weight SNPs.

  - varID\_converter\_fun:
    
    a user defined function to convert varIDs in weights to the
    reference variant format. If NULL, do not convert varIDs.

  - rsID\_converter\_fun:
    
    a user defined function to convert rsIDs in weights to the reference
    variant format. If NULL, do not convert rsIDs.

  - outputdir:
    
    output directory.

  - outname:
    
    name of the output weight file.
