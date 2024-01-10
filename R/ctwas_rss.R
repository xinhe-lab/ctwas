#' Causal inference for TWAS using summary statistics
#' 
#' @param z_gene A data frame with two columns: "id", "z". giving the z scores for genes. 
#' Optionally, a "type" column can also be supplied; this is for using multiple sets of weights
#' 
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#' 
#' @param ld_exprfs A character vector of .`expr` or `.expr.gz` files. One file for
#' one chromosome, in the order of 1 to 22. Therefore, the length of this vector
#' needs to be 22.  `.expr.gz` file is gzip compressed `.expr` files. `.expr` is
#' a matrix of imputed expression values, row is for each sample, column is for
#' each gene. Its sample order is same as in files provided by `.pgenfs`. We also
#' assume corresponding `.exprvar` files are present in the same directory.
#' `.exprvar` files are just tab delimited text files, with columns:
#' \describe{
#'   \item{chrom}{chromosome number, numeric}
#'   \item{p0}{gene boundary position, the smaller value}
#'   \item{p1}{gene boundary position, the larger value}
#'   \item{id}{gene id}
#' }
#' Its rows should be in the same order as the columns for corresponding `.expr`
#' files.
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
#'
#' @param thin The proportion of SNPs to be used for the parameter estimation and initial fine 
#' mapping steps. Smaller \code{thin} parameters reduce runtime at the expense of accuracy. The fine mapping step is rerun using full SNPs 
#' for regions with strong gene signals; see \code{rerun_gene_PIP}.
#'
#' @param prob_single Blocks with probability greater than \code{prob_single} of having 1 or fewer effects will be
#' used for parameter estimation
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
#' "inv_gamma" places an inverse-gamma prior on the variance parameters for each group, with shape and rate hypeparameters.
#' "shared_type" allows all groups in one molecular QTL type to share the same variance parameter.
#' 
#' @param inv_gamma_shape the shape hyperparameter if using "inv_gamma" for \code{group_prior_var_structure}
#' 
#' @param inv_gamma_rate the rate hyperparameter if using "inv_gamma" for \code{group_prior_var_structure}
#'
#' @param estimate_group_prior TRUE/FALSE. If TRUE, the prior inclusion probabilities for SNPs and genes are estimated
#' using the data. If FALSE, \code{group_prior} must be specified
#' 
#' @param estimate_group_prior_var TRUE/FALSE. If TRUE, the prior variances for SNPs and genes are estimated
#' using the data. If FALSE, \code{group_prior_var} must be specified
#' 
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#' 
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#' 
#' @param standardize TRUE/FALSE. If TRUE, all variables are standardized to unit variance
#' 
#' @param ncore The number of cores used to parallelize susie over regions
#' 
#' @param ncore.rerun integer, number of cores to rerun regions with strong signals
#' using full SNPs.
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
#' @param fine_map TRUE/FALSE. If FALSE, only the parameter estimation sets will be run. This is useful for dividing up large jobs
#' into parameter estimation and fine mapping steps
#' 
#' @param reuse_regionlist TRUE/FALSE. If FALSE, call index_regions function to get gene and SNP index and correlation matrix for each region.
#' If True, skip the index_regions function and load pre-computed indexed regions. 
#' @param compress_LDR TRUE/FALSE. If FALSE, correlation matrix folder are not compressed. If TRUE, compressed. 
#' @importFrom logging addHandler loginfo
#' @importFrom tools file_ext
#'
#' @export
#' 
ctwas_rss <- function(
  z_gene,
  z_snp,
  ld_exprvarfs,
  ld_pgenfs = NULL,
  ld_R_dir = NULL,
  ld_regions = c("EUR", "ASN", "AFR"),
  ld_regions_version = c("b37", "b38"),
  ld_regions_custom = NULL,
  thin = 1,
  prob_single = 0.8,
  max_snp_region = Inf,
  rerun_gene_PIP = 0.5,
  niter1 = 3,
  niter2 = 30,
  L= 5,
  group_prior = NULL,
  group_prior_var = NULL,
  group_prior_var_structure = c("independent","shared_all","shared+snps","inv_gamma","shared_type"),
  inv_gamma_shape=1,
  inv_gamma_rate=0,
  estimate_group_prior = T,
  estimate_group_prior_var = T,
  use_null_weight = T,
  coverage = 0.95,
  stardardize = T,
  ncore = 1,
  ncore.rerun = 1,
  ncore_LDR = 1,
  outputdir = getwd(),
  outname = NULL,
  logfile = NULL,
  merge = TRUE,
  fine_map = T,
  reuse_regionlist = F,
  compress_LDR = F){

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
  z_snp$context <- "SNP"   
  if (is.null(z_gene$type)){
    z_gene$type <- "gene"
  }
  if (is.null(z_gene$context)){
    z_gene$context <- "gene"
  }
  
  zdf <- rbind(z_snp[, c("id", "z", "type", "context")], z_gene[, c("id", "z", "type", "context")]) 
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

  if (isTRUE(estimate_group_prior) | isTRUE(estimate_group_prior_var)){

    loginfo("Run susie iteratively, getting rough estimate ...")

    if (!is.null(group_prior)){
      group_prior["SNP"] <- group_prior["SNP"]/thin
    }

    pars <- susieI_rss(zdf = zdf,
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
                       estimate_group_prior = estimate_group_prior,
                       estimate_group_prior_var = estimate_group_prior_var,
                       use_null_weight = use_null_weight,
                       coverage = coverage,
                       ncore = ncore,
                       outputdir = outputdir,
                       outname = paste0(outname, ".s1"),
                       inv_gamma_shape=inv_gamma_shape,
                       inv_gamma_rate=inv_gamma_rate
                   )

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

    pars <- susieI_rss(zdf = zdf,
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
                       estimate_group_prior = estimate_group_prior,
                       estimate_group_prior_var = estimate_group_prior_var,
                       use_null_weight = use_null_weight,
                       coverage = coverage,
                       ncore = ncore,
                       outputdir = outputdir,
                       outname = paste0(outname, ".s2"),
                       inv_gamma_shape=inv_gamma_shape,
                       inv_gamma_rate=inv_gamma_rate)

    group_prior <- pars[["group_prior"]]
    group_prior_var <- pars[["group_prior_var"]]
  }
  
  if (fine_map){
    loginfo("Run susie for all regions.")
    
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
                       estimate_group_prior = estimate_group_prior,
                       estimate_group_prior_var = estimate_group_prior_var,
                       group_prior_var_structure = group_prior_var_structure,
                       use_null_weight = use_null_weight,
                       coverage = coverage,
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
      regionlist <- index_regions(regionfile = regionfile,
                                  exprvarfs = ld_exprvarfs,
                                  pvarfs = ld_pvarfs,
                                  ld_Rfs = ld_Rfs,
                                  select = zdf,
                                  thin = 1, maxSNP = max_snp_region, minvar = 2,
                                  outname = outname, outputdir = outputdir,
                                  merge = merge,
                                  ncore = ncore_LDR,
                                  reuse_R_gene = T) # susie_rss can't take 1 var.
      
      res <- data.table::fread(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"))
      
      # filter out regions based on max gene PIP of the region
      res.keep <- NULL
      for (b in 1: length(regionlist)){
        for (rn in names(regionlist[[b]])){
          #gene_PIP <- max(res$susie_pip[res$type != "SNP" & res$region_tag1 == b & res$region_tag2 == rn], 0)
          gene_PIP <- sum(res$susie_pip[res$type != "SNP" & res$region_tag1 == b & res$region_tag2 == rn])
          if (gene_PIP < rerun_gene_PIP) {
            regionlist[[b]][[rn]] <- NULL
            res.keep <- rbind(res.keep, res[res$region_tag1 ==b & res$region_tag2 == rn, ])
          }
        }
      }
      
      nreg <- sum(unlist(lapply(regionlist, length)))
      
      loginfo("Number of regions that contain strong gene signals: %s", nreg)
      if (nreg == 0){
        file.rename(paste0(file.path(outputdir, outname), ".temp.susieIrss.txt"),
                    paste0(file.path(outputdir, outname), ".susieIrss.txt"))
        
      } else {
        
        loginfo("Rerun susie for regions with strong gene signals using full SNPs.")
        
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
                           estimate_group_prior = estimate_group_prior,
                           estimate_group_prior_var = estimate_group_prior_var,
                           use_null_weight = use_null_weight,
                           coverage = coverage,
                           ncore = ncore.rerun,
                           outputdir = outputdir,
                           outname = paste0(outname, ".s3"),
                           inv_gamma_shape=inv_gamma_shape,
                           inv_gamma_rate=inv_gamma_rate,
                           report_parameters=F)
        
        res.rerun <- data.table::fread(paste0(file.path(outputdir, outname), ".s3.susieIrss.txt"))
        
        res <- rbind(res.keep, res.rerun)
        
        data.table::fwrite(res, file = paste0(file.path(outputdir, outname), ".susieIrss.txt"),
                           sep = "\t", quote = F)
        file.remove((paste0(file.path(outputdir, outname), ".temp.susieIrss.txt")))
      }
    }
  }

  list("group_prior" = group_prior,
       "group_prior_var" = group_prior_var)
  
  if(compress_LDR){
    system(paste0("tar -zcvf ", outputdir, outname, "_LDR.tar.gz ", outputdir, outname, "_LDR"))
  }
}


