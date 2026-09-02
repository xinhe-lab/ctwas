# Loads precomputed correlation matrices for a single region.

It loads precomputed correlation matrices for a single region. It could
load correlation matrices by `region_id` and directory of correlation
matrices `cor_dir`. Otherwise, it loads correlation matrices by the
filenames (`R_sg_file`, `R_sg_file`, `R_s_file`) if they are provided.

## Usage

``` r
load_region_cor(region_id, cor_dir, R_sg_file, R_g_file, R_s_file)
```

## Arguments

  - region\_id:
    
    a character string of region id.

  - cor\_dir:
    
    a string, the directory to store correlation matrices.

  - R\_sg\_file:
    
    filename of SNP-gene correlations.

  - R\_g\_file:
    
    filename of gene-gene correlations.

  - R\_s\_file:
    
    filename of SNP-SNP correlations.

## Value

correlation matrices (R\_snp, R\_snp\_gene and R\_gene)
