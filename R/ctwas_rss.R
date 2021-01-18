#' Causal inference for TWAS using summary statistics
#' @param zdf A data frame with two columns: "id", "z. giving the z scores for
#' genes and snps.
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
#' @param outname a string, the output name
#'
#'
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
ctwas_rss <- function(zdf,
                  ld_pgenfs,
                  ld_exprfs,
                  ld_regions = c("EUR", "ASN", "AFR"),
                  ld_regions_custom = NULL,
                  thin = 1,
                  prob_single = 0.8,
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

  regionlist <- index_regions(ld_pvarfs, ld_exprvarfs, regionfile,
                              select = zdf$id,
                              thin = thin, minvar = 2) # susie_rss can't take 1 var.

  if (isTRUE(estimate_group_prior) | isTRUE(estimate_group_prior_var)){

    loginfo("Run susie iteratively, getting rough estimate ...")

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
                 regionlist = regionlist2,
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

  list("group_prior" = group_prior,
       "group_prior_var" = group_prior_var)

}


