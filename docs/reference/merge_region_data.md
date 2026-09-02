# Merges region data for cross-boundary genes

Merges region data for cross-boundary genes

## Usage

``` r
merge_region_data(
  boundary_genes,
  region_data,
  region_info,
  LD_map,
  snp_map,
  z_snp,
  z_gene,
  expand = TRUE,
  maxSNP = Inf,
  ncore = 1,
  verbose = FALSE,
  logfile = NULL
)
```

## Arguments

  - boundary\_genes:
    
    a data frame of boundary gene info

  - region\_data:
    
    a list of original region\_data

  - region\_info:
    
    a data frame of region definitions

  - LD\_map:
    
    a data frame with filenames of LD matrices and SNP information for
    the regions.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - z\_gene:
    
    A data frame with columns: "id", "z", giving the z-scores for genes.

  - expand:
    
    If TRUE, expand merged region\_data with full SNPs when "thin" \< 1.

  - maxSNP:
    
    Inf or integer. Maximum number of SNPs in a region. Default is Inf,
    no limit. This can be useful if there are many SNPs in a region and
    you don't have enough memory to run the program.

  - ncore:
    
    The number of cores used to parallelize susie over regions

  - verbose:
    
    If TRUE, print detail messages

  - logfile:
    
    The log filename. If NULL, will print log info on screen.

## Value

a list of merged region data, merged region info, LD\_map, snp\_map, and
merged region IDs.
