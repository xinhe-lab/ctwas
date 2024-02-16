#' Screen regions
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param param a list of estimated parameters, including \code{group_prior} and \code{group_prior_var}.
#'
#' @param region_info a data frame of region definition and associated file names
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
#' @param reuse_regionlist TRUE/FALSE. If FALSE, call index_regions function to get gene and SNP index and correlation matrix for each region.
#' If True, skip the index_regions function and load pre-computed indexed regions.
#'
#' @param reuse_regionlist_allSNPs TRUE/FALSE. If FALSE, call index_regions function to get gene and all SNP index and correlation matrix for each region.
#' If True, skip the index_regions function and load pre-computed indexed regions.
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
  param,
  region_info,
  thin = 1,
  max_snp_region = Inf,
  rerun_gene_PIP = 0.5,
  L = 5,
  use_null_weight = T,
  coverage = 0.95,
  min_abs_corr = 0.5,
  ncore = 1,
  outputdir = getwd(),
  outname = NULL,
  logfile = NULL,
  reuse_regionlist = F,
  reuse_regionlist_allSNPs = F,
  compress_LDR = T){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('Screen regions ... ')

  if (length(ld_exprvarfs) != 22){
    stop("Not all imputed expression files for 22 chromosomes are provided.")
  }

  if (is.null(ld_pgenfs) & is.null(ld_R_dir)){
    stop("Error: need to provide either .pgen file or ld_R file")
  } else if (!is.null(ld_pgenfs)){
    if (length(ld_pgenfs) != 22){
      stop("Not all pgen files for 22 chromosomes are provided.")
    }
    ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
    ld_snpinfo <- lapply(ld_pvarfs, read_pvar)
    ld_Rfs <- NULL # do not use R matrix info if genotype is given
  } else {
    ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)
    ld_snpinfo <- lapply(ld_Rfs, read_ld_Rvar)
    ld_pvarfs <- NULL
  }

  if (is.null(ld_regions_custom)){
    ld_regions <- match.arg(ld_regions)
    ld_regions_version <- match.arg(ld_regions_version)
    regionfile <- system.file("extdata", "ldetect",
                              paste0(ld_regions, "." , ld_regions_version, ".bed"), package="ctwas")
  } else {
    regionfile <- ld_regions_custom
  }

  loginfo("LD region file: %s", regionfile)

  z_snp$type <- "SNP"
  z_snp$QTLtype <- "SNP"
  if (is.null(z_gene$type)){
    z_gene$type <- "gene"
  }
  if (is.null(z_gene$QTLtype)){
    z_gene$QTLtype <- "gene"
  }

  zdf <- rbind(z_snp[, c("id", "z", "type", "QTLtype")], z_gene[, c("id", "z", "type", "QTLtype")])
  group_prior_var_structure <- match.arg(group_prior_var_structure)

  rm(z_snp, ld_snpinfo)

  if (thin <=0 | thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  if (!reuse_regionlist){
    regionlist <- index_regions(regionfile = regionfile,
                              exprvarfs = ld_exprvarfs,
                              pvarfs = ld_pvarfs,
                              ld_Rfs = ld_Rfs,
                              select = zdf$id,
                              thin = thin,
                              minvar = 2, # susie_rss can't take 1 var.
                              outname = outname,
                              outputdir = outputdir,
                              merge = merge,
                              ncore = ncore_LDR)

    saveRDS(regionlist, file=paste0(outputdir, "/", outname, ".regionlist.RDS"))

    temp_regs <- lapply(1:22, function(x) cbind(x,
                                              unlist(lapply(regionlist[[x]], "[[", "start")),
                                              unlist(lapply(regionlist[[x]], "[[", "stop"))))

    regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))

    write.table(regs , file= paste0(outputdir,"/", outname, ".regions.txt")
               , row.names=F, col.names=T, sep="\t", quote = F)
  }
  else{
    regionlist <- readRDS(paste0(outputdir, "/", outname, ".regionlist.RDS"))
    regs <- fread(paste0(outputdir,"/", outname, ".regions.txt"))
  }

  loginfo("Run susie for all regions.")

  if (!is.null(group_prior) & !is.null(group_prior_var)){
    load(paste0(outputdir, "/", outname, ".s2.susieIrssres.Rd"))
    group_prior <- group_prior_rec[,ncol(group_prior_rec)]
    group_prior_var <- group_prior_var_rec[,ncol(group_prior_var_rec)]
  }

  pars <- susieI_rss(zdf = zdf,
                     regionlist = regionlist,
                     ld_exprvarfs = ld_exprvarfs,
                     ld_pgenfs = ld_pgenfs,
                     ld_Rfs = ld_Rfs,
                     niter = 1,
                     L = L,
                     z_ld_weight = 0,
                     group_prior = group_prior,
                     group_prior_var = group_prior_var,
                     estimate_group_prior = F,
                     estimate_group_prior_var = F,
                     group_prior_var_structure = group_prior_var_structure,
                     use_null_weight = use_null_weight,
                     coverage = coverage,
                     min_abs_corr = min_abs_corr,
                     ncore = ncore,
                     outputdir = outputdir,
                     outname = paste0(outname, ".temp"),
                     inv_gamma_shape=inv_gamma_shape,
                     inv_gamma_rate=inv_gamma_rate,
                     report_parameters=F)

  group_prior["SNP"] <- group_prior["SNP"] * thin # convert snp pi1

  if (thin == 1) {
    file.rename(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"),
                paste0(file.path(outputdir, outname), ".susieIrss.txt"))
  } else {
    # get full SNPs
    if(!reuse_regionlist_allSNPs){
      regionlist <- index_regions(regionfile = regionfile,
                                  exprvarfs = ld_exprvarfs,
                                  pvarfs = ld_pvarfs,
                                  ld_Rfs = ld_Rfs,
                                  select = zdf,
                                  thin = 1,
                                  maxSNP = max_snp_region,
                                  minvar = 2,
                                  outname = paste0(outname,"_allSNPs"),
                                  outputdir = outputdir,
                                  merge = merge,
                                  ncore = ncore_LDR,
                                  reuse_R_gene = T) # susie_rss can't take 1 var.

      saveRDS(regionlist, file=paste0(outputdir, "/", outname, ".allSNPs.regionlist.RDS"))
      temp_regs <- lapply(1:22, function(x) cbind(x,
                                            unlist(lapply(regionlist[[x]], "[[", "start")),
                                            unlist(lapply(regionlist[[x]], "[[", "stop"))))

      regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
      write.table(regs , file= paste0(outputdir,"/", outname, ".allSNPs.regions.txt")
              , row.names=F, col.names=T, sep="\t", quote = F)
    }
    else{
      regionlist <- readRDS(paste0(outputdir, "/", outname, ".allSNPs.regionlist.RDS"))
      regs <- fread(paste0(outputdir,"/", outname, ".allSNPs.regions.txt"))
    }

    res <- data.table::fread(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"))
    # filter out regions based on max gene PIP of the region
    res.keep <- NULL
    for (b in 1: length(regionlist)){
      for (rn in names(regionlist[[b]])){
        #gene_PIP <- max(res$susie_pip[res$type != "SNP" & res$region_tag1 == b & res$region_tag2 == rn], 0)
        gene_PIP <- sum(res$susie_pip[res$type != "SNP" & res$region_tag1 == b & res$region_tag2 == rn])
        gene_PIP <- ifelse(is.na(gene_PIP), 0, gene_PIP) #0 if gene_PIP is NA (no genes in this region)
        if (gene_PIP < rerun_gene_PIP) {
          regionlist[[b]][[rn]] <- NULL
          res.keep <- rbind(res.keep, res[res$region_tag1 ==b & res$region_tag2 == rn, ])
        }
      }
    }

    nreg <- sum(unlist(lapply(regionlist, length)))

    loginfo("Number of regions that contain strong gene signals: %s", nreg)

  }


  return(regionlist)
}


