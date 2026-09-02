# Loads weights in FUSION format

Loads weights in FUSION format

## Usage

``` r
load_fusion_weights(
  weight_dir,
  fusion_method = c("lasso", "enet", "top1", "blup", "bslmm", "best.cv"),
  fusion_genome_version = NA,
  make_extra_table = FALSE,
  ncore = 1
)
```

## Arguments

  - weight\_dir:
    
    the directory containing FUSION weights ('.wgt.RDat' files).

  - fusion\_method:
    
    a string, specifying the method to choose in FUSION models.
    "best.cv" option will use the best model (smallest p-value) under
    cross-validation.

  - fusion\_genome\_version:
    
    a string, specifying the genome version of FUSION models

  - make\_extra\_table:
    
    If TRUE, make an extra table in predictDB format

  - ncore:
    
    integer, number of cores for parallel computing.
