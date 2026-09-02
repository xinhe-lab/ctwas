# Adds or updates z-scores in region\_data based on z\_snp and z\_gene. this will also update sid and gid based on z\_snp and z\_gene.

Adds or updates z-scores in region\_data based on z\_snp and z\_gene.
this will also update sid and gid based on z\_snp and z\_gene.

## Usage

``` r
update_region_z(
  region_data,
  z_snp,
  z_gene,
  update = c("all", "snps", "genes"),
  ncore = 1
)
```

## Arguments

  - region\_data:
    
    a list of assembled region data.

  - z\_snp:
    
    A data frame with columns: "id", "z", giving the z-scores for SNPs.

  - z\_gene:
    
    A data frame with columns: "id", "z", giving the z-scores for genes.

  - update:
    
    options to update z-scores in region data. "all": update all data
    (default), "snps": updates SNP data only, "genes": updates gene data
    only.

  - ncore:
    
    The number of cores used to parallelize susie over regions
