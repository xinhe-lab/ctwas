# Assembles data for all the regions

Assembles data for all the regions

## Usage

``` r
assemble_region_data(
  region_info,
  z_snp,
  z_gene,
  weights,
  snp_map,
  thin = 1,
  maxSNP = Inf,
  min_group_size = 100,
  trim_by = c("z", "random"),
  thin_by = c("ref", "gwas"),
  adjust_boundary_genes = TRUE,
  ncore = 1,
  seed = 99,
  logfile = NULL
)
```

## Arguments

  - region\_info:
    
    a data frame of region definitions.

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - z\_gene:
    
    A data frame with columns: "id", "z", giving the z-scores for genes.

  - weights:
    
    a list of preprocessed weights.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - thin:
    
    The proportion of SNPs to be used for the parameter estimation and
    initial screening region steps. Smaller `thin` parameters reduce
    runtime at the expense of accuracy. The fine mapping step is rerun
    using full SNPs for regions with strong gene signals.

  - maxSNP:
    
    Inf or integer. Maximum number of SNPs in a region. Default is Inf,
    no limit. This can be useful if there are many SNPs in a region and
    you don't have enough memory to run the program.

  - min\_group\_size:
    
    Minimum number of genes for a group to be included.

  - trim\_by:
    
    remove SNPs if the total number of SNPs exceeds `maxSNP`, options:
    "z" (trim SNPs with lower |z|), "random".

  - thin\_by:
    
    options for thinning SNPs, "reference": thin reference SNPs, "gwas":
    thin GWAS SNPs.

  - adjust\_boundary\_genes:
    
    If TRUE, identify boundary genes, and adjust region\_data for
    boundary genes.

  - ncore:
    
    The number of cores used to parallelize susie over regions

  - seed:
    
    seed for random sampling

  - logfile:
    
    path to the log file, if NULL will print log info on screen.

## Value

a list of region\_data
