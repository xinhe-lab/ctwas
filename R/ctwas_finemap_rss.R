#' cTWAS finemapping using summary statistics
#'
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes.
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param ld_exprvarfs A character vector of `.exprvar` files.
#' `.exprvar` files are tab delimited text files, containing position information for each imputed gene,
#' with columns:
#' \describe{
#'   \item{chrom}{chromosome number, numeric}
#'   \item{p0}{gene boundary position, the smaller value}
#'   \item{p1}{gene boundary position, the larger value}
#'   \item{id}{gene id}
#' }
#'
#' @param ld_pgenfs A character vector of .pgen or .bed files. One file for one
#' chromosome, in the order of 1 to 22. Therefore, the length of this vector
#' needs to be 22. If .pgen files are given, then .pvar and .psam are assumed
#' to present in the same directory. If .bed files are given, then .bim and
#' .fam files are assumed to present in the same directory.
#'
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref".
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#'
#' @param ld_regions A string representing the population to use for defining
#' LD regions. These LD regions were previously defined by ldetect. The user can also
#' provide custom LD regions matching genotype data, see
#' \code{ld_regions_custom}.
#'
#' @param ld_regions_version A string representing the genome reference build ("b37", "b38") to use for defining
#' LD regions. See \code{ld_regions}.
#'
#' @param ld_regions_custom A bed format file defining LD regions. The default
#' is \code{NULL}; when specified, \code{ld_regions} and \code{ld_regions_version} will be ignored.
#'#'
#' @param reuse_regionlist TRUE/FALSE. If FALSE, call index_regions function to get gene and SNP index and correlation matrix for each region.
#' If TRUE, skip the index_regions function and load pre-computed indexed regions.
#'
#' @param regionlist_custom a list object of pre-computed indexed regions.
#'
#' @param region_tags a character vector of region tags, in the format of "chromosome"_"region number". If reuse_regionlist is TRUE,
#' will select a subset of regions from the pre-computed regionlist with the specified region tags.
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
#' "inv_gamma" places an inverse-gamma prior on the variance parameters for each group, with shape and rate hypeparameters.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#'
#' @param inv_gamma_shape the shape hyperparameter if using "inv_gamma" for \code{group_prior_var_structure}
#'
#' @param inv_gamma_rate the rate hyperparameter if using "inv_gamma" for \code{group_prior_var_structure}
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
#' @param ncore_LDR integer, number of cores to use to parallelize construction of the SNP x gene and gene x gene LD matrices
#'
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @param merge TRUE/FALSE. If TRUE, merge regions when a gene spans a region boundary (i.e. belongs to multiple regions.)
#'
#' @param compress_LDR TRUE/FALSE. If FALSE, correlation matrix folder are not compressed. If TRUE, compressed.
#'
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
#'
ctwas_finemap_rss <- function(
    z_gene,
    z_snp,
    ld_exprvarfs,
    ld_pgenfs = NULL,
    ld_R_dir = NULL,
    ld_regions = c("EUR", "ASN", "AFR"),
    ld_regions_version = c("b37", "b38"),
    ld_regions_custom = NULL,
    reuse_regionlist = F,
    regionlist_custom = NULL,
    region_tags = NULL,
    L= 5,
    group_prior = NULL,
    group_prior_var = NULL,
    group_prior_var_structure = c("independent","shared_all","shared+snps","inv_gamma","shared_QTLtype"),
    inv_gamma_shape=1,
    inv_gamma_rate=0,
    use_null_weight = T,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ncore = 1,
    ncore_LDR = 1,
    outputdir = getwd(),
    outname = NULL,
    logfile = NULL,
    merge = TRUE,
    compress_LDR = F){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo('ctwas finemapping regions ... ')

  if (is.null(ld_pgenfs) & is.null(ld_R_dir)){
    stop("Error: need to provide either .pgen file or ld_R file")
  } else if (!is.null(ld_pgenfs)){
    ld_pvarfs <- sapply(ld_pgenfs, prep_pvar, outputdir = outputdir)
    ld_snpinfo <- lapply(ld_pvarfs, read_pvar)
    ld_Rfs <- NULL # do not use R matrix info if genotype is given
  } else {
    ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)
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

  group_prior_var_structure <- match.arg(group_prior_var_structure)

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
                                thin = thin, minvar = 2,
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
        loginfo("select %d regions from the regionlist", length(region_tags))
        subset_regionlist_res <- subset_regionlist(regionlist, region_tags)
        regionlist <- subset_regionlist_res$regionlist
        regs <- subset_regionlist_res$regs
      }

    }
  }

  loginfo("Run finemapping with L = %d", L)

  if (is.null(group_prior) || is.null(group_prior_var)) {
    stop("Error: need to provide both group_prior and group_prior_var")
  }

  if (is.null(regionlist)) {
    stop("Error: need to provide regionlist")
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
                     group_prior_var_structure = group_prior_var_structure,
                     use_null_weight = use_null_weight,
                     coverage = coverage,
                     min_abs_corr = min_abs_corr,
                     ncore = ncore,
                     outputdir = outputdir,
                     outname = outname,
                     inv_gamma_shape=inv_gamma_shape,
                     inv_gamma_rate=inv_gamma_rate,
                     report_parameters=F)

  if(compress_LDR){
    system(paste0("tar -zcvf ", outputdir, outname, "_LDR.tar.gz ", outputdir, outname, "_LDR"))
  }

}

