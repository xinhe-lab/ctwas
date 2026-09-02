# Filter z\_gene by group size

Filter z\_gene by group size

## Usage

``` r
filter_z_gene_by_group_size(z_gene, min_group_size)
```

## Arguments

  - z\_gene:
    
    a data frame of gene z-scores, with columns: "id", "z", "type",
    "context", "group".

  - min\_group\_size:
    
    Minimum number of variables in a group.

## Value

a data frame of gene z-scores.
