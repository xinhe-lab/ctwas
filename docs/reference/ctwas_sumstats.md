# cTWAS analysis using summary statistics

cTWAS analysis using summary statistics

## Usage

``` r
ctwas_sumstats(
  z_snp,
  weights,
  region_info,
  LD_map,
  snp_map,
  z_gene = NULL,
  thin = 1,
  niter_prefit = 3,
  niter = 50,
  L = 5,
  init_group_prior = NULL,
  init_group_prior_var = NULL,
  group_prior_var_structure = c("shared_all", "shared_type", "shared_context",
    "shared_nonSNP", "independent"),
  min_nonSNP_PIP = 0.5,
  min_snp_pval = 5e-08,
  min_p_single_effect = 0.8,
  maxSNP = Inf,
  min_var = 2,
  min_gene = 1,
  min_group_size = 100,
  null_method = c("ctwas", "susie", "none"),
  EM_tol = 1e-04,
  coverage = 0.95,
  min_abs_corr = 0.1,
  include_prior = FALSE,
  include_susie_result = FALSE,
  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
  LD_loader_fun = NULL,
  snpinfo_loader_fun = NULL,
  force_compute_cor = FALSE,
  save_cor = FALSE,
  cor_dir = NULL,
  outputdir = NULL,
  outname = "ctwas",
  ncore = 1,
  ncore_LD = max(ncore - 1, 1),
  seed = 99,
  logfile = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

  - z\_snp:
    
    A data frame with four columns: "id", "A1", "A2", "z". giving the z
    scores for SNPs. "A1" is effect allele. "A2" is the other allele.

  - weights:
    
    a list of pre-processed prediction weights.

  - region\_info:
    
    a data frame of region definitions.

  - LD\_map:
    
    a data frame with filenames of LD matrices and SNP information for
    the regions.

  - snp\_map:
    
    a list of data frames with SNP-to-region map for the reference.

  - z\_gene:
    
    A data frame with columns: "id", "z", giving the z-scores for genes.

  - thin:
    
    The proportion of SNPs to be used for estimating parameters and
    screening regions.

  - niter\_prefit:
    
    the number of iterations of the E-M algorithm to perform during the
    initial parameter estimation step.

  - niter:
    
    the maximum number of iterations of the E-M algorithm to perform
    during the complete parameter estimation step.

  - L:
    
    the number of effects for susie during the fine mapping steps.

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

  - min\_nonSNP\_PIP:
    
    Regions with non-SNP PIP \>= `min_nonSNP_PIP` will be selected to
    run finemapping using all SNPs.

  - min\_snp\_pval:
    
    Select regions with minimum GWAS SNP p-values \< `min_snp_pval`.

  - min\_p\_single\_effect:
    
    Regions with probability greater than `min_p_single_effect` of
    having 1 or fewer effects will be used for parameter estimation.

  - maxSNP:
    
    Inf or integer. Maximum number of SNPs in a region. Default is Inf,
    no limit. This can be useful if there are many SNPs in a region and
    you don't have enough memory to run the program.

  - min\_var:
    
    minimum number of variables (SNPs and genes) in a region when
    estimating paramters and screening regions.

  - min\_gene:
    
    minimum number of genes in a region when estimating paramters and
    screening regions.

  - min\_group\_size:
    
    Minimum number of genes in a group. Groups with number of genes \<
    `min_group_size` will be removed for the analysis.

  - null\_method:
    
    Method to compute null model, options: "ctwas", "susie" or "none".

  - EM\_tol:
    
    A small, non-negative number specifying the convergence tolerance of
    log-likelihood for the EM iterations.

  - coverage:
    
    A number between 0 and 1 specifying the “coverage” of the estimated
    confidence sets.

  - min\_abs\_corr:
    
    Minimum absolute correlation allowed in a credible set.

  - include\_prior:
    
    If TRUE, include priors in finemapping results.

  - include\_susie\_result:
    
    If TRUE, include the "susie" result object in finemapping results.

  - LD\_format:
    
    file format for LD matrix. If "custom", use a user defined
    `LD_loader_fun()` function to load LD matrix.

  - LD\_loader\_fun:
    
    a user defined function to load LD matrix when `LD_format =
    "custom"`.

  - snpinfo\_loader\_fun:
    
    a user defined function to load SNP information file, if SNP
    information files are not in standard cTWAS reference format.

  - force\_compute\_cor:
    
    If TRUE, force computing correlation (R) matrices.

  - save\_cor:
    
    If TRUE, save correlation (R) matrices to `cor_dir`.

  - cor\_dir:
    
    The directory to store correlation (R) matrices.

  - outputdir:
    
    The directory to store output. If specified, save outputs to the
    directory.

  - outname:
    
    The output name.

  - ncore:
    
    The number of cores used to parallelize computing over regions.

  - ncore\_LD:
    
    The number of cores used to parallelize computing correlation
    matrices, in screening regions and fine-mapping steps with LD.

  - seed:
    
    seed for random sampling when thinning the SNPs in region data.

  - logfile:
    
    The log filename. If NULL, print log info on screen.

  - verbose:
    
    If TRUE, print detailed messages.

  - ...:
    
    Additional arguments of `susie_rss`.

## Value

a list, including z\_gene, estimated parameters, region\_data, screening
region results, and fine-mapping results.
