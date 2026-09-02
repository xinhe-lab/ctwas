# Converts PLINK genotype data to LD matrices and SNP info files, saves LD matrices as .RDS files and SNP info as .Rvar files

Converts PLINK genotype data to LD matrices and SNP info files, saves LD
matrices as .RDS files and SNP info as .Rvar files

## Usage

``` r
convert_geno_to_LD_matrix(
  region_info,
  genotype_files,
  varinfo_files,
  chrom = 1:22,
  outputdir = getwd(),
  outname = "",
  include_variance = TRUE,
  include_allele_freq = TRUE,
  show_progress_bar = TRUE,
  verbose = FALSE,
  logfile = NULL
)
```

## Arguments

  - region\_info:
    
    a data frame of region definitions, with columns: chrom, start,
    stop, and region\_id.

  - genotype\_files:
    
    Reference genotype files in PLINK binary genotype data in .pgen or
    .bed format. It should contain files for all chromosomes (from 1 to
    22), one file per chromosome.

  - varinfo\_files:
    
    Reference variant information files in PLINK .pvar or .bim format.
    It could have one file per chromosome or have one big file for all
    chromosomes. The output will use the genome positions in
    `varinfo_files`.

  - chrom:
    
    a vector of chromosome numbers to process genotype data.

  - outputdir:
    
    Output directory.

  - outname:
    
    Output filestem.

  - include\_variance:
    
    If TRUE, include variance in .Rvar output.

  - include\_allele\_freq:
    
    If TRUE, include allele frequency in .Rvar output.

  - show\_progress\_bar:
    
    If TRUE, print progress bar.

  - verbose:
    
    If TRUE, print detail messages.

  - logfile:
    
    The log filename. If NULL, print log info on screen.

## Value

a data frame of region\_metatable, with region definitions and filenames
of LD matrices and variant information.
