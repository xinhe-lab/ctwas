#' Screen regions
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weight_list
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param regionlist a list object indexing regions, variants and genes.
#'
#' @param weight_list a list of weights
#'
#' @param thin The proportion of SNPs to be used for parameter estimation and initial screening regions.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs
#' for regions with strong gene signals; see \code{min_gene_PIP}.
#'
#' @param max_snp_region Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the last rerun step
#' (using full SNPs and rerun susie for regions with strong gene signals) only.
#'
#' @param min_gene_PIP Regions with non-SNP PIP >= \code{min_gene_PIP}
#' will be selected to run finemapping using full SNPs.
#'
#' @param L the number of effects for susie
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of regionlist, screened region tags, fine mapping result for regions with weak signal
#'
#' @export
#'
screen_regions <- function(
    z_snp,
    z_gene,
    region_info,
    gene_info = NULL,
    weight_list = NULL,
    regionlist = NULL,
    thin = 1,
    max_snp_region = Inf,
    min_gene_PIP = 0.5,
    L = 1,
    group_prior = NULL,
    group_prior_var = NULL,
    use_null_weight = TRUE,
    coverage = 0.95,
    min_abs_corr = 0.5,
    max_iter = 1,
    ncore = 1,
    logfile = NULL,
    ...){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  # combine z-scores of SNPs and genes
  z_snp$type <- "SNP"
  z_snp$QTLtype <- "SNP"
  if (is.null(z_gene$type)){
    z_gene$type <- "gene"
  }
  if (is.null(z_gene$QTLtype)){
    z_gene$QTLtype <- "gene"
  }
  zdf <- rbind(z_snp[, c("id", "z", "type", "QTLtype")],
               z_gene[, c("id", "z", "type", "QTLtype")])

  if (thin <= 0 || thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  loginfo("Screen regions ...")

  if (is.null(regionlist)) {
    loginfo("Get regionlist with thin = %.2f", thin)
    res <- get_regionlist(region_info = region_info,
                          gene_info = gene_info,
                          weight_list = weight_list,
                          select = zdf$id,
                          thin = thin,
                          maxSNP = max_snp_region,
                          minvar = 2)

    regionlist <- res$regionlist
    rm(res)
  }

  # run finemapping for all regions containing thinned SNPs
  loginfo("Run initial screening for all regions ...")

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  corelist <- region2core(regionlist, ncore)

  finemap_res <- foreach (core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {
    susie_res.core.list <- list()
    # run susie for each region
    region_tags.core <- corelist[[core]]
    for (region_tag in region_tags.core) {
      susie_res <- finemap_region(z_snp = z_snp,
                                  z_gene = z_gene,
                                  gene_info = gene_info,
                                  regionlist = regionlist,
                                  region_tag = region_tag,
                                  weight_list = weight_list,
                                  L = L,
                                  group_prior = group_prior,
                                  group_prior_var = group_prior_var,
                                  use_null_weight = use_null_weight,
                                  coverage = coverage,
                                  min_abs_corr = min_abs_corr,
                                  max_iter = max_iter,
                                  ...)
      susie_res.core.list[[i]] <- susie_res
    }
    susie_res.core <- do.call(rbind, susie_res.core.list)
    susie_res.core
  }
  parallel::stopCluster(cl)

  # select regions based on max gene PIP of the region
  # TODO: confirm the change to use TOTAL none SNP PIP rather than MAX gene PIP?
  finemap_weak_res <- NULL
  screened_region_tags <- NULL
  screened_region_info <- NULL
  for (region_tag in names(regionlist)){
    #gene_PIP <- max(finemap_res$susie_pip[finemap_res$type != "SNP" & finemap_res$region_tag == region_tag], 0)
    region_finemap_res <- finemap_res[finemap_res$region_tag == region_tag,]
    gene_PIP <- sum(region_finemap_res$susie_pip[region_finemap_res$type != "SNP"])
    gene_PIP[is.na(gene_PIP)] <- 0 # 0 if gene_PIP is NA (no genes in this region)
    if (gene_PIP >= min_gene_PIP) {
      screened_region_tags <- c(screened_region_tags, region_tag)
      screened_region_info <- rbind(screened_region_info, region_info[region_info$region_tag = region_tag, ])
    }
  }

  loginfo("Number of region tags that contain strong gene signals: %d", length(screened_region_tags))

  # get regionlist containing all SNPs for screened regions
  if(thin == 1){
    screened_regionlist <- regionlist[screened_region_tags]
  }else{
    # TODO: update regionlist with screened region tags
    screened_regionlist <- update_regionlist(...)
  }

  # keep the finemapping results for the regions without strong signals (will not rerun)
  weak_region_finemap_res <- finemap_res[!finemap_res$region_tag %in% screened_region_tags, ]

  return(list("screened_regionlist" = updated_regionlist,
              "screened_region_tags" = screened_region_tags,
              "weak_region_finemap_res" = weak_region_finemap_res))
}


