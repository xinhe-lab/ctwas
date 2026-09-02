# Runs EM to estimate parameters using cTWAS SER model

Runs EM to estimate parameters using cTWAS SER model

## Usage

``` r
fit_EM(
  region_data,
  groups,
  types,
  contexts,
  niter = 30,
  init_group_prior = NULL,
  init_group_prior_var = NULL,
  group_prior_var_structure = c("shared_all", "shared_type", "shared_context",
    "shared_nonSNP", "independent", "fixed"),
  null_method = c("ctwas", "susie", "none"),
  EM_tol = 1e-04,
  force_run_niter = FALSE,
  warn_converge_fail = TRUE,
  ncore = 1,
  verbose = FALSE
)
```

## Arguments

  - region\_data:
    
    a list object with the susie input data for each region

  - groups:
    
    a vector of the groups to estimate parameters

  - types:
    
    a vector of the types to estimate parameters

  - contexts:
    
    a vector of the contexts to estimate parameters

  - niter:
    
    the maximum number of iterations of the EM algorithm to perform

  - init\_group\_prior:
    
    a vector of initial prior inclusion probabilities for SNPs and
    genes.

  - init\_group\_prior\_var:
    
    a vector of initial prior variances for SNPs and gene effects.

  - group\_prior\_var\_structure:
    
    a string indicating the structure to put on the prior variance
    parameters. "shared\_all" allows all groups to share the same
    variance parameter. "shared\_type" allows all groups in one
    molecular QTL type to share the same variance parameter.
    "shared\_context" allows all groups in one context (tissue, cell
    type, condition) to share the same variance parameter.
    "shared\_nonSNP" allows all non-SNP groups to share the same
    variance parameter. "independent" allows all groups to have their
    own separate variance parameters.

  - null\_method:
    
    Method to compute null weight, options: "ctwas", "susie" or "none".

  - EM\_tol:
    
    A small, non-negative number specifying the convergence tolerance of
    log-likelihood for the EM iterations.

  - force\_run\_niter:
    
    If TRUE, run all the `niter` iterations.

  - warn\_converge\_fail:
    
    If `warn_converge_fail = TRUE`, prints a warning message when EM
    algorithm does not converge.

  - ncore:
    
    The number of cores used to parallelize over regions

  - verbose:
    
    If TRUE, print detail messages

## Value

a list of estimated parameters
