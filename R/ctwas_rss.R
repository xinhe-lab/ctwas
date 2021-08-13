#' Causal inference for TWAS using summary statistics
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#' If `harmonize= False`, A1 and A2 are not required.
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for
#' genes.
#'
#' @param ld_pgenfs A character vector of .pgen or .bed files. One file for one
#'  chromosome, in the order of 1 to 22. Therefore, the length of this vector
#'  needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#'  to present in the same directory. If .bed files are given, then .bim and
#'  .fam files are assumed to present in the same directory.
#'
#' @param ld_exprfs A character vector of .`expr` or `.expr.gz` files. One file for
#'  one chromosome, in the order of 1 to 22. Therefore, the length of this vector
#'  needs to be 22.  `.expr.gz` file is gzip compressed `.expr` files. `.expr` is
#'  a matrix of imputed expression values, row is for each sample, column is for
#'  each gene. Its sample order is same as in files provided by `.pgenfs`. We also
#'  assume corresponding `.exprvar` files are present in the same directory.
#'  `.exprvar` files are just tab delimited text files, with columns:
#'  \describe{
#'    \item{chrom}{chromosome number, numeric}
#'    \item{p0}{gene boundary position, the smaller value}
#'    \item{p1}{gene boundary position, the larger value}
#'    \item{id}{gene id}
#'  }
#'  Its rows should be in the same order as the columns for corresponding `.expr`
#'  files.
#'
#' @param ld_regions A string representing the population to use for defining
#'  LD regions. These LD regions are defined by ldetect. The user can also
#'  provide custom LD regions matching genotype data, see
#'  \code{ld_regions_custom}.
#'
#' @param ld_regions_custom A bed format file defining LD regions. The default
#'  is \code{NULL}; when specified, \code{ld_regions} will be ignored.
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param prob_single blocks with probility more than this value will be
#'  used for parameter estimation
#'
#' @param rerun_gene_PIP if thin <1, will rerun blocks with the max gene PIP
#' > \code{rerun_gene_PIP} using full SNPs. if \code{rerun_gene_PIP} is 0, then
#' all blocks will rerun with full SNPs
#'
#' @param harmonize TRUE/FALSE. If TRUE, will harmonize GWAS, LD reference
#' and weight internally (flip alleles and remove strand ambiguous SNPs)
#'
#' @param max_snp_region Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program. This applies to the last rerun step
#'  (using full SNPs and rerun susie for regions with strong gene signals) only.
#'
#' @param ncore integer, number of cores to run parameter estimation
#' @param ncore.rerun integer, number of cores to rerun regions with strong signals
#' using full SNPs.
#'
#' @param outname a string, the output name
#'
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
ctwas_rss <- function(
  z_gene,
  z_snp,
  ld_exprfs,
  ld_pgenfs = NULL,
  ld_R_dir = NULL,
  ld_regions = c("EUR", "ASN", "AFR"),
  ld_regions_version = c("b37", "b38"),
  ld_regions_custom = NULL,
  thin = 1,
  prob_single = 0.8,
  rerun_gene_PIP = 0.8,
  niter1 = 3,
  niter2 = 30,
  L= 5,
  group_prior = NULL,
  group_prior_var = NULL,
  estimate_group_prior = T,
  estimate_group_prior_var = T,
  use_null_weight = T,
  coverage = 0.95,
  stardardize = T,
  max_snp_region = Inf,
  ncore = 1,
  ncore.rerun = 1,
  outputdir = getwd(),
  outname = NULL,
  logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('ctwas started ... ')

  if (length(ld_exprfs) != 22){
    stop("Not all imputed expression files for 22 chromosomes are provided.")
  }
  ld_exprvarfs <- sapply(ld_exprfs, prep_exprvar)

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

  zdf <- rbind(z_snp[, c("id", "z")], z_gene[, c("id", "z")])
  rm(z_snp, ld_snpinfo)

  if (thin <=0 | thin > 1){
    stop("thin value needs to be in (0,1]")
  }

  regionlist <- index_regions(regionfile = regionfile,
                              exprvarfs = ld_exprvarfs,
                              pvarfs = ld_pvarfs,
                              ld_Rfs = ld_Rfs,
                              select = zdf$id,
                              thin = thin, minvar = 2,
                              outname = outname,
                              outputdir = outputdir) # susie_rss can't take 1 var.
  
  temp_regs <- lapply(1:22, function(x) cbind(x,
                   unlist(lapply(regionlist[[x]], "[[", "start")),
                     unlist(lapply(regionlist[[x]], "[[", "stop"))))
                 
  regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))

  write.table(regs , file= paste0(outputdir,"/", outname, ".regions.txt")
               , row.names=F, col.names=T, sep="\t", quote = F)

  if (isTRUE(estimate_group_prior) | isTRUE(estimate_group_prior_var)){

    loginfo("Run susie iteratively, getting rough estimate ...")

    if (!is.null(group_prior)){
      group_prior[2] <- group_prior[2]/thin
    }

    pars <- susieI_rss(zdf = zdf,
                       regionlist = regionlist,
                       ld_exprfs = ld_exprfs,
                       ld_pgenfs = ld_pgenfs,
                       ld_Rfs = ld_Rfs,
                       niter = niter1,
                       L = 1,
                       z_ld_weight = 0,
                       group_prior = group_prior,
                       group_prior_var = group_prior_var,
                       estimate_group_prior = estimate_group_prior,
                       estimate_group_prior_var = estimate_group_prior_var,
                       use_null_weight = use_null_weight,
                       coverage = coverage,
                       ncore = ncore,
                       outputdir = outputdir,
                       outname = paste0(outname, ".s1")
                   )


    group_prior <- pars[["group_prior"]]
    group_prior_var <- pars[["group_prior_var"]]

    # filter blocks
    regionlist2 <- filter_regions(regionlist,
                                  group_prior,
                                  prob_single= prob_single)

    loginfo("Blocks are filtered: %s blocks left",
            sum(unlist(lapply(regionlist2, length))))

    loginfo("Run susie iteratively, getting accurate estimate ...")

    pars <- susieI_rss(zdf = zdf,
                       regionlist = regionlist2,
                       ld_exprfs = ld_exprfs,
                       ld_pgenfs = ld_pgenfs,
                       ld_Rfs = ld_Rfs,
                       niter = niter2,
                       L = 1,
                       z_ld_weight = 0,
                       group_prior = group_prior,
                       group_prior_var = group_prior_var,
                       estimate_group_prior = estimate_group_prior,
                       estimate_group_prior_var = estimate_group_prior_var,
                       use_null_weight = use_null_weight,
                       coverage = coverage,
                       ncore = ncore,
                       outputdir = outputdir,
                       outname = paste0(outname, ".s2"))


    group_prior <- pars[["group_prior"]]
    group_prior_var <- pars[["group_prior_var"]]
  }

  loginfo("Run susie for all regions.")

  pars <- susieI_rss(zdf = zdf,
                     regionlist = regionlist,
                     ld_exprfs = ld_exprfs,
                     ld_pgenfs = ld_pgenfs,
                     ld_Rfs = ld_Rfs,
                     niter = 1,
                     L = L,
                     z_ld_weight = 0,
                     group_prior = group_prior,
                     group_prior_var = group_prior_var,
                     estimate_group_prior = estimate_group_prior,
                     estimate_group_prior_var = estimate_group_prior_var,
                     use_null_weight = use_null_weight,
                     coverage = coverage,
                     ncore = ncore,
                     outputdir = outputdir,
                     outname = paste0(outname, ".temp"))


  group_prior[2] <- group_prior[2] * thin # convert snp pi1

  if (thin == 1) {
    file.rename(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"),
                paste0(file.path(outputdir, outname), ".susieIrss.txt"))

  } else {

    # get full SNPs
    regionlist <- index_regions(regionfile = regionfile,
                                exprvarfs = ld_exprvarfs,
                                pvarfs = ld_pvarfs,
                                ld_Rfs = ld_Rfs,
                                select = zdf,
                                thin = 1, maxSNP = max_snp_region, minvar = 2,
                                outname = outname, outputdir = outputdir) # susie_rss can't take 1 var.

    res <- data.table::fread(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"))

    # filter out regions based on max gene PIP of the region
    res.keep <- NULL
    for (b in 1: length(regionlist)){
      for (rn in names(regionlist[[b]])){
       gene_PIP <- max(res[res$type == "gene" & res$region_tag1 == b & res$region_tag2 == rn, ]$susie_pip, 0)
       if (gene_PIP < rerun_gene_PIP) {
         regionlist[[b]][[rn]] <- NULL
         res.keep <- rbind(res.keep, res[res$region_tag1 ==b & res$region_tag2 == rn, ])
       }
      }
    }

    nreg <- sum(unlist(lapply(regionlist, length)))

    loginfo("Number of regions that contains strong gene signals: %s", nreg)
    if (nreg == 0){
      file.rename(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"),
                  paste0(file.path(outputdir, outname), ".susieIrss.txt"))

    } else {

      loginfo("Rerun susie for regions with strong gene signals using full SNPs.")

      pars <- susieI_rss(zdf = zdf,
                         regionlist = regionlist,
                         ld_exprfs = ld_exprfs,
                         ld_pgenfs = ld_pgenfs,
                         ld_Rfs = ld_Rfs,
                         niter = 1,
                         L = L,
                         z_ld_weight = 0,
                         group_prior = group_prior,
                         group_prior_var = group_prior_var,
                         estimate_group_prior = estimate_group_prior,
                         estimate_group_prior_var = estimate_group_prior_var,
                         use_null_weight = use_null_weight,
                         coverage = coverage,
                         ncore = ncore.rerun,
                         outputdir = outputdir,
                         outname = paste0(outname, ".s3"))

      res.rerun <- data.table::fread(paste0(file.path(outputdir, outname), ".s3.susieIrss.txt"))

      res <- rbind(res.keep, res.rerun)

      data.table::fwrite(res, file = paste0(file.path(outputdir, outname), ".susieIrss.txt"),
                         sep = "\t", quote = F)
      file.remove((paste0(file.path(outputdir, outname), ".temp.susieIrss.txt")))
    }
  }

  list("group_prior" = group_prior,
       "group_prior_var" = group_prior_var)

}


