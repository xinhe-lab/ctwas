#' cTWAS fine-map specific regions with full SNPs (thin = 1)
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param regionlist the list of regions to be fine-mapped.
#'
#' @param region_tags a character string of region tags
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
#' @return a list of finemapping results.
#'
#' @export
#'
finemap_regions <- function(z_gene,
                             z_snp,
                             regionlist,
                             region_info,
                             region_tags = NULL,
                             L = 5,
                             group_prior = NULL,
                             group_prior_var = NULL,
                             use_null_weight = T,
                             coverage = 0.95,
                             min_abs_corr = 0.5,
                             ncore = 1,
                             logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  loginfo("Finemapping regions with all SNPs... ")

  # combine z-scores of different types
  zdf <- combine_z(z_gene, z_snp)

  # Compute correlation matrices and generate regionlist
  if (is.null(regionlist)){
    loginfo("Compute correlation matrices and generate regionlist with all SNPs (thin = 1)")

    regionlist <- compute_cor(region_info,
                              gene_info,
                              weight_list,
                              select = zdf$id,
                              thin = 1,
                              minvar = 2,
                              outname = outname,
                              outputdir = outputdir,
                              merge = merge,
                              ncore = ncore)

    saveRDS(regionlist, file=paste0(outputdir, "/", outname, ".regionlist.RDS"))

    # temp_regs <- lapply(1:22, function(x) cbind(x,
    #                                             unlist(lapply(regionlist[[x]], "[[", "start")),
    #                                             unlist(lapply(regionlist[[x]], "[[", "stop"))))
    #
    # regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
    #
    # write.table(regs , file= paste0(outputdir,"/", outname, ".regions.txt")
    #             , row.names=F, col.names=T, sep="\t", quote = F)
  } else {
    loginfo("Use existing regionlist")
  }

  # select and assemble a subset of regionlist by region_tags
  if (length(region_tags) > 0){
    loginfo("subset %s regions", length(region_tags))
    subset_regionlist_res <- subset_regionlist(regionlist, region_tags)
    regionlist <- subset_regionlist_res$regionlist
  }

  loginfo("Run cTWAS finemapping with L = %d", L)

  # run finemapping using SuSiE RSS
  finemap_res <- ctwas_susieI_rss(zdf = zdf,
                                  regionlist = regionlist,
                                  region_info = region_info,
                                  niter = 1,
                                  L = L,
                                  group_prior = group_prior,
                                  group_prior_var = group_prior_var,
                                  estimate_group_prior = FALSE,
                                  estimate_group_prior_var = FALSE,
                                  use_null_weight = use_null_weight,
                                  coverage = coverage,
                                  min_abs_corr = min_abs_corr,
                                  ncore = ncore,
                                  verbose = F)

  return(finemap_res)

}

