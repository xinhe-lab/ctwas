# Make convergence plots for the estimated parameters

Make convergence plots for the estimated parameters

## Usage

``` r
make_convergence_plots(
  param,
  gwas_n,
  plot_loglik = FALSE,
  log_enrichment = FALSE,
  title.size = 10,
  legend.size = 8,
  colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
    "#999999"),
  group_names,
  ncol = 2
)
```

## Arguments

  - param:
    
    a list of cTWAS parameter estimation result from `est_param`

  - gwas\_n:
    
    the sample size of the GWAS summary statistics

  - plot\_loglik:
    
    If TRUE, plot log-likelihood.

  - log\_enrichment:
    
    If TRUE, plot enrichment in log scale.

  - title.size:
    
    font size of the plot title

  - legend.size:
    
    font size of the plot legend title

  - colors:
    
    colors for different groups, passed as the `values` input to
    `scale_color_manual`. If fewer colors than "fits" are given, the
    colors are recycled.

  - group\_names:
    
    groups to be included in the plots

  - ncol:
    
    number of columns in the output plot
