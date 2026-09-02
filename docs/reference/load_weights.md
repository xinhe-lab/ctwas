# Loads weights in PredictDB or FUSION format

Loads weights in PredictDB or FUSION format

## Usage

``` r
load_weights(
  weight_file,
  weight_format = c("PredictDB", "FUSION"),
  filter_protein_coding_genes = TRUE,
  load_predictdb_LD = TRUE,
  fusion_method = c("lasso", "enet", "top1", "blup", "bslmm", "best.cv"),
  fusion_genome_version = NA,
  ncore = 1
)
```

## Arguments

  - weight\_file:
    
    filename of the '.db' file for PredictDB weights; or the directory
    containing '.wgt.RDat' files for FUSION weights.

  - weight\_format:
    
    a string or a vector, specifying format of each weight file, e.g.
    PredictDB, FUSION.

  - filter\_protein\_coding\_genes:
    
    If TRUE, keep protein coding genes only. This option is only for
    PredictDB weights

  - load\_predictdb\_LD:
    
    If TRUE, load pre-computed LD among weight SNPs. This option is only
    for PredictDB weights

  - fusion\_method:
    
    a string, specifying the method to choose in FUSION models.
    "best.cv" option will use the best model (smallest p-value) under
    cross-validation.

  - fusion\_genome\_version:
    
    a string, specifying the genome version of FUSION models

  - ncore:
    
    integer, number of cores for parallel computing.
