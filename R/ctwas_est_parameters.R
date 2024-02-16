#' Estimate cTWAS parameters
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param thin The proportion of SNPs to be used for the parameter estimation and initial fine
#' mapping steps. Smaller \code{thin} parameters reduce runtime at the expense of accuracy.
#' The fine mapping step is rerun using full SNPs for regions with strong gene signals.
#'
#' @param prob_single Blocks with probability greater than \code{prob_single} of having 1 or fewer effects will be
#' used for parameter estimation
#'
#' @param niter1 the number of iterations of the E-M algorithm to perform during the initial parameter estimation step
#'
#' @param niter2 the number of iterations of the E-M algorithm to perform during the complete parameter estimation step
#'
#' @param L the number of effects for susie during the fine mapping steps
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes. This is ignored
#' if \code{estimate_group_prior = T}
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects. This is ignored
#' if \code{estimate_group_prior_var = T}
#'
#' @param group_prior_var_structure a string indicating the structure to put on the prior variance parameters.
#' "independent" is the default and allows all groups to have their own separate variance parameters.
#' "shared" allows all groups to share the same variance parameter.
#' "shared+snps" allows all groups to share the same variance parameter, and this variance parameter is also shared with SNPs.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
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
#' @return a list with estimated parameters, and updated region_info: containing correlation file names for each region.
#'
#' @export
#'
est_param <- function(
    z_gene,
    z_snp,
    region_info,
    thin = 1,
    prob_single = 0.8,
    niter1 = 3,
    niter2 = 30,
    L= 5,
    group_prior = NULL,
    group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared+snps","shared_QTLtype"),
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

  loginfo('ctwas started ... ')

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
                                thin = thin, minvar = 2,
                                outname = outname,
                                outputdir = outputdir,
                                merge = merge,
                                ncore = ncore_LDR) # susie_rss can't take 1 var.

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

  loginfo("Run susie iteratively, getting rough estimate ...")

  pars <- ctwas_susieI_rss(zdf = zdf,
                           regionlist = regionlist,
                           ld_exprvarfs = ld_exprvarfs,
                           ld_pgenfs = ld_pgenfs,
                           ld_Rfs = ld_Rfs,
                           niter = niter1,
                           L = 1,
                           z_ld_weight = 0,
                           group_prior = group_prior,
                           group_prior_var = group_prior_var,
                           group_prior_var_structure = group_prior_var_structure,
                           estimate_group_prior = T,
                           estimate_group_prior_var = T,
                           use_null_weight = use_null_weight,
                           coverage = coverage,
                           min_abs_corr = min_abs_corr,
                           ncore = ncore,
                           outputdir = outputdir,
                           outname = paste0(outname, ".s1"))

  group_prior <- pars[["group_prior"]]
  group_prior_var <- pars[["group_prior_var"]]

  # filter blocks
  regionlist2 <- filter_regions(regionlist,
                                group_prior,
                                prob_single,
                                zdf)

  loginfo("Blocks are filtered: %s blocks left",
          sum(unlist(lapply(regionlist2, length))))

  loginfo("Run susie iteratively, getting accurate estimate ...")

  pars <- ctwas_susieI_rss(zdf = zdf,
                           regionlist = regionlist2,
                           ld_exprvarfs = ld_exprvarfs,
                           ld_pgenfs = ld_pgenfs,
                           ld_Rfs = ld_Rfs,
                           niter = niter2,
                           L = 1,
                           z_ld_weight = 0,
                           group_prior = group_prior,
                           group_prior_var = group_prior_var,
                           group_prior_var_structure = group_prior_var_structure,
                           estimate_group_prior = T,
                           estimate_group_prior_var = T,
                           use_null_weight = use_null_weight,
                           coverage = coverage,
                           min_abs_corr = min_abs_corr,
                           ncore = ncore,
                           outputdir = outputdir,
                           outname = paste0(outname, ".s2"))
}

