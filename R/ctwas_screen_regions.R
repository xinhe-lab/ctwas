#' Screen regions
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param weights a string, pointing to a directory with the FUSION/TWAS format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param thin The proportion of SNPs to be used for parameter estimation and initial screening regions.
#' Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs
#' for regions with strong gene signals; see \code{rerun_gene_PIP}.
#'
#' @param max_snp_region Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the last rerun step
#' (using full SNPs and rerun susie for regions with strong gene signals) only.
#'
#' @param rerun_gene_PIP if thin <1, will rerun blocks with the max gene PIP
#' > \code{rerun_gene_PIP} using full SNPs. if \code{rerun_gene_PIP} is 0, then
#' all blocks will rerun with full SNPs
#'
#' @param L the number of effects for susie during the fine mapping steps
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
    weights = NULL,
    regionlist = NULL,
    regionlist_allSNPs = NULL,
    thin = 1,
    max_snp_region = Inf,
    rerun_gene_PIP = 0.5,
    L = 1,
    group_prior = NULL,
    group_prior_var = NULL,
    use_null_weight = TRUE,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ncore = 1,
    logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  # combine z-scores of different types
  zdf <- combine_z(z_snp, z_gene)

  if (thin <= 0 || thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  # get regionlist containing all SNPs (thin = 1)
  if (is.null(regionlist_allSNPs)) {
    loginfo("Get regionlist containing all SNPs")
    res <- get_region_idx(region_info = region_info,
                          gene_info = gene_info,
                          weights = weights,
                          select = zdf$id,
                          thin = 1,
                          maxSNP = max_snp_region,
                          minvar = 2,
                          merge = FALSE,
                          ncore = ncore)

    regionlist_full <- res$regionlist
    rm(res)
  }

  # TODO: confirm if thin = 1, skip screening regions
  if (thin != 1) {

    loginfo("thin = 1, skip screening regions.")
    regionlist <- regionlist_allSNPs
    screened_region_tags <- get_region_tags(regionlist_allSNPs)
    finemap_weak_res <- NULL

  }else{

    loginfo("Screen regions ...")

    # get regionlist containing thinned SNPs
    if (is.null(regionlist)) {
      loginfo("Get regionlist with thin = %.2f", thin)
      res <- get_region_idx(region_info = region_info,
                            gene_info = gene_info,
                            weights = weights,
                            select = zdf$id,
                            thin = thin,
                            minvar = 2,
                            merge = FALSE,
                            ncore = ncore)

      regionlist_thinned <- res$regionlist
      rm(res)
    }

    # run finemapping for all regions containing thinned SNPs
    loginfo("Run initial screening for all regions ...")

    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)

    corelist <- region2core(regionlist_thinned, ncore)

    finemap_res <- foreach (core = 1:length(corelist), .combine = "rbind", .packages = "ctwas") %dopar% {

      susie_res.core.list <- list()
      # run susie for each region
      regs <- corelist[[core]]
      for (i in 1: nrow(regs)) {
        b <- regs[i, "b"]
        rn <- regs[i, "rn"]
        temp_region_tag <- paste0(b, ":", rn)
        susie_res <- finemap_region(z_snp = z_snp,
                                    z_gene = z_gene,
                                    region_info = region_info,
                                    gene_info = gene_info,
                                    regionlist = regionlist_thinned,
                                    region_tag = temp_region_tag,
                                    L = L,
                                    group_prior = group_prior,
                                    group_prior_var = group_prior_var,
                                    use_null_weight = use_null_weight,
                                    coverage = coverage,
                                    min_abs_corr = min_abs_corr)
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
    # screened_regionlist <- regionlist_allSNPs
    for (b in 1: length(regionlist_allSNPs)){
      for (rn in names(regionlist_allSNPs[[b]])){
        #gene_PIP <- max(finemap_res$susie_pip[finemap_res$type != "SNP" & finemap_res$region_tag1 == b & finemap_res$region_tag2 == rn], 0)
        gene_PIP <- sum(finemap_res$susie_pip[finemap_res$type != "SNP" & finemap_res$region_tag1 == b & finemap_res$region_tag2 == rn])
        gene_PIP <- ifelse(is.na(gene_PIP), 0, gene_PIP) # 0 if gene_PIP is NA (no genes in this region)
        if (gene_PIP < rerun_gene_PIP) {
          # screened_regionlist[[b]][[rn]] <- NULL
          # keep the finemapping results for the regions without strong signals (will not rerun)
          finemap_weak_res <- rbind(finemap_weak_res, finemap_res[finemap_res$region_tag1 == b & finemap_res$region_tag2 == rn, ])
        }else{
          region_tag <- paste0(b, ":", rn)
          screened_region_tags <- c(screened_region_tags, region_tag)
        }
      }
    }

    # nreg <- sum(unlist(lapply(screened_regionlist, length)))
    # loginfo("Number of regions that contain strong gene signals: %d", nreg)
    loginfo("Number of regions that contain strong gene signals: %d", length(screened_region_tags))

  }

  return(list("region_info" = region_info,
              "regionlist" = regionlist_thinned,
              "regionlist_allSNPs" = regionlist_allSNPs,
              "screened_region_tags" = screened_region_tags,
              "finemap_weak_res" = finemap_weak_res))
}


