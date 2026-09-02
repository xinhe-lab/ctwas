# Estimates cTWAS parameters using EM with cTWAS SER model

Estimates cTWAS parameters using EM with cTWAS SER model

## Usage

``` r
est_param(
  region_data,
  init_group_prior = NULL,
  init_group_prior_var = NULL,
  group_prior_var_structure = c("shared_all", "shared_type", "shared_context",
    "shared_nonSNP", "independent", "fixed"),
  niter_prefit = 3,
  niter = 50,
  min_var = 2,
  min_gene = 1,
  min_group_size = 100,
  min_p_single_effect = 0.8,
  null_method = c("ctwas", "susie", "none"),
  EM_tol = 1e-04,
  force_run_niter = FALSE,
  ncore = 1,
  logfile = NULL,
  verbose = FALSE
)
```

## Arguments

  - region\_data:
    
    a list object indexing regions, variants and genes.

  - init\_group\_prior:
    
    a vector of initial values of prior inclusion probabilities for
    different groups.

  - init\_group\_prior\_var:
    
    a vector of initial values of prior variances for different groups.

  - group\_prior\_var\_structure:
    
    a string indicating the structure to put on the prior variance
    parameters. "shared\_all" allows all groups to share the same
    variance parameter. "shared\_type" allows all groups in one
    molecular QTL type to share the same variance parameter.
    "shared\_context" allows all groups in one context (tissue, cell
    type, condition) to share the same variance parameter.
    "shared\_nonSNP" allows all non-SNP groups to share the same
    variance parameter. "independent" allows all groups to have their
    own separate variance parameters. "fixed" sets prior variance
    parameters to values in `init_group_prior_var`.

  - niter\_prefit:
    
    the maximum number of iterations of the E-M algorithm to perform
    during the initial parameter estimation step.

  - niter:
    
    the maximum number of iterations of the E-M algorithm to perform
    during the complete parameter estimation step.

  - min\_var:
    
    minimum number of variables (SNPs and genes) in a region.

  - min\_gene:
    
    minimum number of genes in a region.

  - min\_group\_size:
    
    Minimum number of variables in a group.

  - min\_p\_single\_effect:
    
    Regions with probability greater than `min_p_single_effect` of
    having at most one causal effect will be used selected for the
    complete parameter estimation step.

  - null\_method:
    
    Method to compute null model, options: "ctwas", "susie" or "none".

  - EM\_tol:
    
    A small, non-negative number specifying the convergence tolerance of
    log-likelihood for the EM iterations.

  - force\_run\_niter:
    
    If TRUE, run all the `niter` EM iterations.

  - ncore:
    
    The number of cores used to parallelize computation over regions.

  - logfile:
    
    The log filename. If NULL, print log info on screen.

  - verbose:
    
    If TRUE, print detail messages.

## Value

a list with estimated parameters
