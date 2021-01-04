# ctwas
package for the causal TWAS project

# Install

`ctwas` depends on R package `pgenlibr`. We need to install `pgenlibr` first.

```
install.packages("devtools")
devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
```

or 

```
$ git clone https://github.com/chrchang/plink-ng.git
$ R
> install.packages('plink-ng/2.0/pgenlibr', repos = NULL, type='source')
```

Install `ctwas`.

```
devtools::install_github("simingz/ctwas", ref = "main")
```

# Notes

SuSiE R code obtained from:
```
git clone git@github.com:stephenslab/susieR.git
git checkout a06e35e
```

Minor edits to make it accept different prior variances for each variable. 
