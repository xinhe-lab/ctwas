# Preprocess PredictDB/FUSION weights and harmonize with LD reference

Preprocess PredictDB/FUSION weights and harmonize with LD reference

## Usage

``` r
preprocess_weights(
  weight_file,
  region_info,
  gwas_snp_ids,
  snp_map,
  LD_map = NULL,
  type,
  context,
  weight_name = paste0(context, "_", type),
  weight_format = c("PredictDB", "FUSION"),
  drop_strand_ambig = TRUE,
  filter_protein_coding_genes = TRUE,
  scale_predictdb_weights = TRUE,
  load_predictdb_LD = TRUE,
  include_weight_LD = TRUE,
  fusion_method = c("lasso", "enet", "top1", "blup", "bslmm", "best.cv"),
  fusion_genome_version = NA,
  top_n_snps = NULL,
  LD_format = c("rds", "rdata", "csv", "txt", "custom"),
  LD_loader_fun = NULL,
  snpinfo_loader_fun = NULL,
  varID_converter_fun = NULL,
  ncore = 1,
  logfile = NULL,
  verbose = FALSE
)
```

## Arguments

  - weight\_file:
    
    filename of the '.db' file for PredictDB weights; or the directory
    containing '.wgt.RDat' files for FUSION weights.

  - region\_info:
    
    a data frame of region definitions.

  - gwas\_snp\_ids:
    
    a vector of SNP IDs in GWAS summary statistics (z\_snp$id).

  - snp\_map:
    
    a list of SNP-to-region map for the reference.

  - LD\_map:
    
    a data frame with filenames of LD matrices and SNP information for
    the regions. Required when `load_predictdb_LD = FALSE`.

  - type:
    
    a string, specifying QTL type of the weight file, e.g. expression,
    splicing, protein.

  - context:
    
    a string, specifying context (tissue/cell type) of the weight file,
    e.g. Liver, Lung, Brain.

  - weight\_name:
    
    a string, specifying name of the weight file. By default, it is
    `weight_name = paste0(context, "_", type)`

  - weight\_format:
    
    a string, specifying format of the weight file, e.g. PredictDB,
    FUSION.

  - drop\_strand\_ambig:
    
    If TRUE remove strand ambiguous variants (A/T, G/C).

  - filter\_protein\_coding\_genes:
    
    If TRUE, keep protein coding genes only. This option is only for
    PredictDB weights.

  - scale\_predictdb\_weights:
    
    If TRUE, scale PredictDB weights by the variance. This is because
    PredictDB weights assume that variant genotypes are not
    standardized, but our implementation assumes standardized variant
    genotypes. This option is only for PredictDB weights.

  - load\_predictdb\_LD:
    
    If TRUE, load pre-computed LD among weight SNPs. This option is only
    for PredictDB weights.

  - include\_weight\_LD:
    
    If TRUE, include LD of variants in weights (R\_wgt) in the weights
    object. R\_wgt is used for computing gene Z-scores. If FALSE, will
    skip computing R\_wgt. This could save running time when using
    precomputed gene Z-scores.

  - fusion\_method:
    
    a string, specifying the method to choose in FUSION models.
    "best.cv" option will use the best model (smallest p-value) under
    cross-validation.

  - fusion\_genome\_version:
    
    a string, specifying the genome version of FUSION models

  - top\_n\_snps:
    
    a number, only keeping the top n SNPs included in weight models. By
    default, keep all SNPs in weights.

  - LD\_format:
    
    file format for LD matrix. If "custom", use a user defined
    `LD_loader_fun()` function to load LD matrix.

  - LD\_loader\_fun:
    
    a user defined function to load LD matrix when `LD_format =
    "custom"`.

  - snpinfo\_loader\_fun:
    
    a user defined function to load SNP information file, if SNP
    information files are not in standard cTWAS reference format.

  - varID\_converter\_fun:
    
    a user defined function to convert weight variant IDs to the
    reference variant format.

  - ncore:
    
    The number of cores used to parallelize computation.

  - logfile:
    
    The log filename. If NULL, will print log info on screen.

  - verbose:
    
    If TRUE, print detail messages.

## Value

a list of processed weights
