# Converts fusion weights to predictDB format

Converts fusion weights to predictDB format

## Usage

``` r
convert_fusion_to_predictdb(
  weight_dir,
  fusion_method = c("lasso", "enet", "top1", "blup", "bslmm", "best.cv"),
  fusion_genome_version = NA,
  make_extra_table = TRUE,
  cov_table = NULL,
  outputdir = getwd(),
  outname
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

  - cov\_table:
    
    a data frame of covariances between variants, with columns:
    "GENE","RSID1","RSID2","VALUE". If NULL, do not create covariance
    files (.txg.gz).

  - outputdir:
    
    output directory

  - outname:
    
    name of the output weight file
