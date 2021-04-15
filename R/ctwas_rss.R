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
#' @param outname a string, the output name
#'
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
ctwas_rss <- function(z_snp,
                  z_gene,
                  ld_pgenfs,
                  ld_exprfs,
                  ld_regions = c("EUR", "ASN", "AFR"),
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
                  harmonize =T,
                  ncore = 1,
                  outputdir = getwd(),
                  outname = NULL,
                  logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('ctwas started ... ')

  if (length(ld_pgenfs) != 22){
    stop("Not all pgen files for 22 chromosomes are provided.")
  }

  if (length(ld_exprfs) != 22){
    stop("Not all imputed expression files for 22 chromosomes are provided.")
  }

  if (is.null(ld_regions_custom)){
    ld_regions <- match.arg(ld_regions)
    regionfile <- system.file("extdata", "ldetect",
                              paste0(ld_regions, ".bed"), package="ctwas")
  } else {
    regionfile <- ld_regions_custom
  }

  loginfo("LD region file: %s", regionfile)

  ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
  ld_exprvarfs <- sapply(ld_exprfs, prep_exprvar)

  if (isTRUE(harmonize)){
    loginfo("flipping z scores to match LD reference")
    for (ld_pvarf in ld_pvarfs){
      ld_snpinfo <- read_pvar(ld_pvarf)
      z_snp <- harmonize_z_ld(z_snp, ld_snpinfo)
    }
  }

  zdf <- rbind(z_snp[, c("id", "z")], z_gene[, c("id", "z")])
  rm(z_snp)

  regionlist <- index_regions(ld_pvarfs, ld_exprvarfs, regionfile,
                              select = zdf$id,
                              thin = thin, minvar = 2) # susie_rss can't take 1 var.

  regs <- do.call(rbind, lapply(1:22, function(x) cbind(x,
                    unlist(lapply(regionlist[[x]], "[[", "start")),
                    unlist(lapply(regionlist[[x]], "[[", "stop")))))

  write.table(regs , file= paste0(outputdir,"/", outname, ".regions.txt")
               , row.names=F, col.names=T, sep="\t", quote = F)

  if (isTRUE(estimate_group_prior) | isTRUE(estimate_group_prior_var)){

    loginfo("Run susie iteratively, getting rough estimate ...")

    if (!is.null(group_prior)){
      group_prior[2] <- group_prior[2]/thin
    }

    pars <- susieI_rss(zdf = zdf,
                       ld_pgenfs = ld_pgenfs,
                       ld_exprfs = ld_exprfs,
                       regionlist = regionlist,
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
                   ld_pgenfs = ld_pgenfs,
                   ld_exprfs = ld_exprfs,
                   regionlist = regionlist2,
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
                 ld_pgenfs = ld_pgenfs,
                 ld_exprfs = ld_exprfs,
                 regionlist = regionlist,
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
                 outname = outname)

  group_prior[2] <- group_prior[2] * thin # convert snp pi1

  if (thin < 1){

    # get full SNPs
    regionlist <- index_regions(ld_pvarfs, ld_exprvarfs, regionfile,
                                select = zdf$id,
                                thin = 1, minvar = 2) # susie_rss can't take 1 var.

    res <- data.table::fread(paste0(file.path(outputdir, outname), ".susieIrss.txt"))

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

    if (nreg > 0){

      loginfo("Rerun susie for regions with strong gene signals using full SNPs.")

      pars <- susieI_rss(zdf = zdf,
                         ld_pgenfs = ld_pgenfs,
                         ld_exprfs = ld_exprfs,
                         regionlist = regionlist,
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
                         outname = paste0(outname, ".s3"))

      res.rerun <- data.table::fread(paste0(file.path(outputdir, outname), ".s3.susieIrss.txt"))

      res <- rbind(res.keep, res.rerun)

      data.table::fwrite(res, file = paste0(file.path(outputdir, outname), ".susieIrss.txt"),
                         sep = "\t", quote = F)
    }
  }

  list("group_prior" = group_prior,
       "group_prior_var" = group_prior_var)

}


