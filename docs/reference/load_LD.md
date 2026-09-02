# Load LD matrix

Load LD matrix

## Usage

``` r
load_LD(
  file,
  format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
  LD_loader_fun = NULL
)
```

## Arguments

  - file:
    
    LD matrix file name.

  - format:
    
    file format for LD matrix. If "custom", use a user defined
    `LD_loader_fun()` function to load LD matrix.

  - LD\_loader\_fun:
    
    a user defined function to load LD matrix
