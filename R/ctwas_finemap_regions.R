#' cTWAS fine-map specific regions with full SNPs (thin = 1)
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
#' @param regionlist the list of regions to be fine-mapped.
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
#'#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom tools file_ext
#'
#' @return a list of PIPs and CS.
#'
#' @export
#'
finemap_regions <- function(
    z_gene,
    z_snp,
    param,
    region_info,
    regionlist,
    L= 5,
    use_null_weight = T,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ncore = 1,
    logfile){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo("cTWAS finemapping regions with full SNPs... ")

  if (is.null(ld_pgenfs) & is.null(ld_R_dir)){
    stop("Error: need to provide either .pgen file or ld_R file")
  } else if (!is.null(ld_pgenfs)){
    ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
    ld_snpinfo <- lapply(ld_pvarfs, read_pvar)
    ld_Rfs <- NULL # do not use R matrix info if genotype is given
  } else {
    ld_Rfs <- file.path(outputdir, paste0(outname, "_ld_R_chr", 1:22, ".txt"))
    if (!all(file.exists(ld_Rfs))) {
      ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)
    }
    ld_snpinfo <- lapply(ld_Rfs, read_ld_Rvar)
    ld_pvarfs <- NULL
  }

  # prepare input z_snp and z_gene data
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
  rm(z_snp, ld_snpinfo)

  if (!reuse_regionlist){
    # index regions using region file and and save regionlist
    if (!is.null(ld_regions_custom)){
      regionfile <- ld_regions_custom
    } else {
      ld_regions <- match.arg(ld_regions)
      ld_regions_version <- match.arg(ld_regions_version)
      regionfile <- system.file("extdata", "ldetect",
                                paste0(ld_regions, "." , ld_regions_version, ".bed"), package="ctwas")
    }
    loginfo("use LD region file: %s", regionfile)

    regionlist <- index_regions(regionfile = regionfile,
                                exprvarfs = ld_exprvarfs,
                                pvarfs = ld_pvarfs,
                                ld_Rfs = ld_Rfs,
                                select = zdf$id,
                                thin = 1,
                                minvar = 2,
                                outname = outname,
                                outputdir = outputdir,
                                merge = merge,
                                ncore = ncore_LDR) # susie_rss can't take 1 var.

    saveRDS(regionlist, file=paste0(outputdir, "/", outname, ".regionlist.RDS"))

    # save regions.txt
    reg <- read.table(regionfile, header = T, stringsAsFactors = F)
    if (is.character(reg$chr)){
      reg$chr <- readr::parse_number(reg$chr)
    }
    reg_chrs <- sort(unique(reg$chr))

    temp_regs <- lapply(reg_chrs, function(x) cbind(x,
                                                    unlist(lapply(regionlist[[x]], "[[", "start")),
                                                    unlist(lapply(regionlist[[x]], "[[", "stop"))))

    regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
    write.table(regs, file= paste0(outputdir,"/", outname, ".regions.txt"),
                row.names=F, col.names=T, sep="\t", quote = F)
  }else{
    if (!is.null(regionlist_custom)) {
      loginfo("use custom regionlist")
      regionlist <- regionlist_custom
    } else {
      loginfo("reuse regionlist: %s", paste0(outputdir, "/", outname, ".regionlist.RDS"))
      regionlist <- readRDS(paste0(outputdir, "/", outname, ".regionlist.RDS"))
      regs <- data.table::fread(paste0(outputdir,"/", outname, ".regions.txt"))
      # select and assemble a subset of regionlist by region_tags
      if (length(region_tags) > 0){
        loginfo("subset %s regions from precomputed regionlist", length(region_tags))
        subset_regionlist_res <- subset_regionlist(regionlist, region_tags)
        regionlist <- subset_regionlist_res$regionlist
        regs <- subset_regionlist_res$regs
      }
    }
  }

  loginfo("Run cTWAS finemapping with L = %d", L)

  if (is.null(group_prior) || is.null(group_prior_var)) {
    stop("Error: need to provide both group_prior and group_prior_var")
  }

  if (is.null(regionlist)) {
    stop("Error: need to provide regionlist")
  }

  # run finemapping using SuSiE RSS
  pars <- ctwas_susieI_rss(zdf = zdf,
                           regionlist = regionlist,
                           ld_exprvarfs = ld_exprvarfs,
                           ld_pgenfs = ld_pgenfs,
                           ld_Rfs = ld_Rfs,
                           niter = 1,
                           L = L,
                           z_ld_weight = 0,
                           group_prior = group_prior,
                           group_prior_var = group_prior_var,
                           use_null_weight = use_null_weight,
                           coverage = coverage,
                           min_abs_corr = min_abs_corr,
                           ncore = ncore,
                           report_parameters=F)

  if(compress_LDR){
    system(paste0("tar -zcvf ", outputdir, outname, "_LDR.tar.gz ", outputdir, outname, "_LDR"))
  }

}

