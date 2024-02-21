#' Screen regions
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param gene_info a data frame of gene information obtained from \code{compute_gene_z}
#'
#' @param weight a string, pointing to a directory with the FUSION/TWAS format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param thin The proportion of SNPs to be used for the parameter estimation and initial fine
#' mapping steps. Smaller \code{thin} parameters reduce runtime at the expense of accuracy. The fine mapping step is rerun using full SNPs
#' for regions with strong gene signals; see \code{rerun_gene_PIP}.

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
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param compress_LDR TRUE/FALSE. If FALSE, correlation matrix folder are not compressed. If TRUE, compressed.
#'
#' @importFrom logging addHandler loginfo
#'
#' @importFrom tools file_ext
#'
#' @return a list of regions to be fine-mapped
#'
#' @export
#'
screen_regions <- function(
  z_gene,
  z_snp,
  region_info,
  gene_info = NULL,
  weight = NULL,
  regionlist = NULL,
  regionlist_allSNPs = NULL,
  thin = 1,
  max_snp_region = Inf,
  rerun_gene_PIP = 0.5,
  group_prior = NULL,
  group_prior_var = NULL,
  L = 5,
  use_null_weight = T,
  coverage = 0.95,
  min_abs_corr = 0.5,
  ncore = 1,
  outputdir = getwd(),
  outname = NULL,
  logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Screen regions for finemapping ... ')

  # combine z-scores of different types
  zdf <- combine_z(z_gene, z_snp)

  if (thin <= 0 || thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  if (thin == 1) {
    loginfo("thin = 1, skip screening regions.")
    return(regionlist)

  } else {
    loginfo("Screening regions.")

    if (is.null(regionlist)){
      loginfo("Compute correlation matrices and generate regionlist with thin = %.2f", thin)

      regionlist <- compute_cor(region_info = region_info,
                                gene_info = gene_info,
                                weight = weight,
                                select = zdf$id,
                                thin = thin,
                                minvar = 2,
                                outname = outname,
                                outputdir = outputdir,
                                merge = FALSE,
                                ncore = ncore)

      saveRDS(regionlist, file=paste0(outputdir, "/", outname, ".regionlist.RDS"))

      # temp_regs <- lapply(1:22, function(x) cbind(x,
      #                                             unlist(lapply(regionlist[[x]], "[[", "start")),
      #                                             unlist(lapply(regionlist[[x]], "[[", "stop"))))
      #
      # regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
      #
      # write.table(regs , file= paste0(outputdir,"/", outname, ".regions.txt"),
      #             row.names=F, col.names=T, sep="\t", quote = F)
    }

    # run finemapping for all regions with thinned correlation matrices
    susieI_res <- ctwas_susieI_rss(zdf = zdf,
                                   region_info = region_info,
                                   regionlist = regionlist,
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
    finemap_res <- susieI_res$susieIrss_res

    group_prior["SNP"] <- group_prior["SNP"] * thin # convert snp pi1

    # get regionlist with all SNPs
    if(is.null(regionlist_allSNPs)){
      loginfo("Compute correlation matrices and generate regionlist with all SNPs (thin = 1)")

      regionlist <- compute_cor(region_info = region_info,
                                gene_info = gene_info,
                                weight = weight,
                                select = zdf$id,
                                thin = 1,
                                maxSNP = max_snp_region,
                                minvar = 2,
                                outname = paste0(outname,".allSNPs"),
                                outputdir = outputdir,
                                merge = FALSE,
                                ncore = ncore,
                                reuse_R_gene = T)

      saveRDS(regionlist, file=paste0(outputdir, "/", outname, ".allSNPs.regionlist.RDS"))
      # temp_regs <- lapply(1:22, function(x) cbind(x,
      #                                             unlist(lapply(regionlist[[x]], "[[", "start")),
      #                                             unlist(lapply(regionlist[[x]], "[[", "stop"))))
      #
      # regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
      # write.table(regs , file= paste0(outputdir,"/", outname, ".allSNPs.regions.txt"),
      #             row.names=F, col.names=T, sep="\t", quote = F)
    }else{
      regionlist <- regionlist_allSNPs
    }

    # filter out regions based on max gene PIP of the region
    res.keep <- NULL
    for (b in 1: length(regionlist)){
      for (rn in names(regionlist[[b]])){
        #gene_PIP <- max(finemap_res$susie_pip[res$type != "SNP" & finemap_res$region_tag1 == b & finemap_res$region_tag2 == rn], 0)
        gene_PIP <- sum(finemap_res$susie_pip[res$type != "SNP" & finemap_res$region_tag1 == b & finemap_res$region_tag2 == rn])
        gene_PIP <- ifelse(is.na(gene_PIP), 0, gene_PIP) #0 if gene_PIP is NA (no genes in this region)
        if (gene_PIP < rerun_gene_PIP) {
          regionlist[[b]][[rn]] <- NULL
          res.keep <- rbind(res.keep, res[res$region_tag1 ==b & res$region_tag2 == rn, ])
        }
      }
    }

    nreg <- sum(unlist(lapply(regionlist, length)))

    loginfo("Number of regions that contain strong gene signals: %d", nreg)

  }

  return(regionlist)
}


