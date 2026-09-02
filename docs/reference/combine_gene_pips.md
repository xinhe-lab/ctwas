# Combines gene PIPs by context, type or group.

Combines gene PIPs by context, type or group.

## Usage

``` r
combine_gene_pips(
  susie_alpha_res,
  mapping_table = NULL,
  map_by = "molecular_id",
  drop_unmapped = TRUE,
  group_by = "molecular_id",
  by = c("context", "type", "group"),
  method = c("combine_cs", "sum"),
  filter_cs = FALSE,
  keep_alpha_in_cs_only = FALSE,
  include_cs_id = TRUE,
  include_set_id = FALSE,
  missing_value = NA
)
```

## Arguments

  - susie\_alpha\_res:
    
    a data frame of annotated susie alpha result.

  - mapping\_table:
    
    a data frame of mapping between molecular traits and genes, with
    required columns: "molecular\_id", "gene\_name".

  - map\_by:
    
    column name to be mapped by (default: "molecular\_id").

  - drop\_unmapped:
    
    If TRUE, remove unmapped genes.

  - group\_by:
    
    column name to group genes by.

  - by:
    
    option to combine PIPs by: "context" (default), "type", or "group".

  - method:
    
    method to combine PIPs of molecular traits targeting the same gene.
    options: "combine\_cs" (default): first sums PIPs of molecular
    traits of a genes in each credible set, and then combine PIPs using
    the following formula: \\(1 - \\prod\_k (1 - \\text{PIP}\_k)\\),
    where \\(\\text{PIP}\_k\\) is the summed PIP of the \\(k\\)-th
    credible set of a gene. This is the default option for combining
    PIPs from fine-mapping with LD. "sum": sum over PIPs of all
    molecular traits for the same gene. This summation is the expected
    number of causal molecular traits in this gene, and could be higher
    than 1. We will use this option for combining PIPs from fine-mapping
    without LD.

  - filter\_cs:
    
    If TRUE, limits gene results to credible sets (CS).

  - keep\_alpha\_in\_cs\_only:
    
    If TRUE, only keep single effects (alpha) in credible sets when
    calculating combined PIP. This is similar to `prune_by_cs` in susie.

  - include\_cs\_id:
    
    If TRUE, include credible set IDs of the genes in the output

  - include\_set\_id:
    
    If TRUE, include susie set IDs of the genes in the output

  - missing\_value:
    
    set missing value as (default: NA)

## Value

a data frame of combined gene PIPs for each context, type or group
