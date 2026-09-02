# Creates weight files in PredictDB format from QTL data

Creates weight files in PredictDB format from QTL data

## Usage

``` r
create_predictdb_from_QTLs(
  weight_table,
  gene_table = NULL,
  cov_table = NULL,
  use_top_QTL = TRUE,
  select_by = c("pval", "weight"),
  outputdir = getwd(),
  outname
)
```

## Arguments

  - weight\_table:
    
    a data frame of the genes, QTLs and weights, with columns: "gene",
    "rsid", "varID", "ref\_allele", "eff\_allele", "weight". If you want
    to use multiple eQTLs per gene, you can set `use_top_QTL=FALSE`. But
    we assume the weights of the eQTLs are learned from multiple
    regression (instead of marginal effect sizes).

  - gene\_table:
    
    a data frame (optional) with information of the genes in
    `weight_table` ("gene","genename","gene\_type", etc.). If NULL,
    create a simple gene\_table based on the weight\_table

  - cov\_table:
    
    a data frame of covariances between variants, with columns:
    "GENE","RSID1","RSID2", "VALUE". If NULL, do not create covariance
    files (.txg.gz), unless `use_top_QTL=TRUE`.

  - use\_top\_QTL:
    
    If TRUE, only keep the top QTL per gene (molecular trait), and
    create a simple cov\_table with covariance set to 1.

  - select\_by:
    
    Select the top SNP by the column: "pval": choose the top SNP with
    the smallest p-value per gene (molecular trait). "weight": choose
    the top SNP with the largest abs(weight) per gene (molecular trait),
    Only used when `use_top_QTL=TRUE`.

  - outputdir:
    
    output directory.

  - outname:
    
    name of the output weight file.
