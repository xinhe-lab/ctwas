#' Causal inference for TWAS using summary statistics
ctwas_sumstats <- function(
    z_snp,
    weight,
    region_info,
    weight_format = c("PredictDB", "FUSION"),
    method = c("lasso", "blup", "bslmm", "top1", "enet", "best"),
    niter1 = 3,
    niter2 = 30,
    L = 5,
    group_prior = NULL,
    group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared+snps","shared_QTLtype"),
    thin = 1,
    prob_single = 0.8,
    use_null_weight = T,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ncore = 1,
    outputdir = getwd(),
    outname = NULL,
    logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  if (thin <= 0 | thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  # compute gene z-scores
  res <- compute_gene_z(z_snp = z_snp,
                        weight = weight,
                        region_info = region_info,
                        weight_format = weight_format,
                        method = method,
                        ncore = ncore)
  z_gene <- res$z_gene
  z_snp <- res$z_snp
  gene_info <- res$gene_info
  rm(res)

  # estimate parameters
  # including computing correlation matrices for thinned SNPs
  res <- est_param(z_snp = z_snp,
                   region_info = region_info,
                   z_gene = z_gene,
                   gene_info = gene_info,
                   thin = thin,
                   group_prior = group_prior,
                   group_prior_var = group_prior_var,
                   group_prior_var_structure = group_prior_var_structure,
                   niter1 = niter1,
                   niter2 = niter2,
                   outputdir = outputdir,
                   outname = outname,
                   ncore = ncore)

  param <- res$param
  group_prior <- param$group_prior
  group_prior_var <- param$group_prior_var
  region_info <- res$region_info
  regionlist <- res$regionlist
  rm(res)

  # screen regions
  # including computing correlation matrices for all SNPs, and finemapping for all regions with thinned SNPs,
  # and selecting regions with strong gene signals
  res <- screen_regions(z_snp = z_snp,
                        z_gene = z_gene,
                        region_info = region_info,
                        gene_info = gene_info,
                        weight = weight,
                        regionlist = regionlist,
                        thin = thin,
                        max_snp_region = max_snp_region,
                        rerun_gene_PIP = rerun_gene_PIP,
                        L = L,
                        group_prior = group_prior,
                        group_prior_var = group_prior_var,
                        use_null_weight = use_null_weight,
                        coverage = coverage,
                        min_abs_corr = min_abs_corr,
                        outputdir = outputdir,
                        outname = outname,
                        ncore = ncore)

  region_info <- res$region_info
  screened_regionlist <- res$screened_regionlist
  finemap_weak_res <- res$finemap_weak_res
  rm(res)

  # convert SNP prior
  group_prior["SNP"] <- group_prior["SNP"] * thin
  param$group_prior <- group_prior

  # run fine-mapping for regions with strong gene signals using all SNPs
  finemap_strong_res <- finemap_regions(z_snp = z_snp,
                                        z_gene = z_gene,
                                        region_info = region_info,
                                        regionlist = screened_regionlist,
                                        L = L,
                                        group_prior = group_prior,
                                        group_prior_var = group_prior_var,
                                        use_null_weight = use_null_weight,
                                        coverage = coverage,
                                        min_abs_corr = min_abs_corr)

  # combine fine-mapping results for regions with weak signals and strong signals
  finemap_res <- rbind(finemap_weak_res, finemap_strong_res)

  return(list("param" = param,
              "region_info" = region_info,
              "z_gene" = z_gene,
              "finemap_res" = finemap_res))

}

