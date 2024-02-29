#' Get gene and SNP index for each region
#' @description For each region, get the index for snp and gene (index
#' is location/column number in .pgen file or .expr file) located within
#' this region.
#'
#' @param regionfile regions file. Has three columns: chr, start, end. The regions file
#' should provide non overlaping regions defining LD blocks. currently does not support
#' chromsome X/Y etc.
#'
#' @param select Default is NULL, all variants will be selected. Or a vector of variant IDs,
#' or a data frame with columns id and z (id is for gene or SNP id, z is for z scores).
#' z will be used for remove SNPs if the total number of SNPs exceeds limit. See
#' parameter `maxSNP` for more information.
#' 
#' @param thin  A scalar in (0,1]. The proportion of SNPs
#' left after down sampling. Only applied on SNPs after selecting variants.
#'  
#' @param maxSNP Default is Inf, no limit for the maximum number of SNPs in a region. Or an
#' integer indicating the maximum number of SNPs allowed in a region. This
#' parameter is useful when a region contains many SNPs and you don't have enough memory to
#' run the program. In this case, you can put a limit on the number of SNPs in the region.
#' If z scores are given in the parameter `select`, i.e. a data frame with columns id and z is
#' provided, SNPs are ranked based on |z| from high to low and only the top `maxSNP` SNPs
#' are kept. If only variant ids are provided, then `maxSNP` number of SNPs will be chosen
#' randomly.
#'  
#' @param minvar minimum number of variatns in a region
#'
#' @param merge TRUE/FALSE. If TRUE, merge regions when a gene spans a region boundary (i.e. belongs to multiple regions.)
#' 
#' @param outputdir a string, the directory to store output
#' 
#' @param outname a string, the output name
#' 
#' @param ncore the number of cores used to parallelize region indexing
#' 
#' @param reuse_R_gene an option to reuse the R_gene matrix when indexing for the final rerun step
#'
#' @return A list. Items correspond to each pvarf/exprvarf. Each Item is
#'  also a list, the items in this list are for each region.
#'
#' @importFrom logging loginfo
#'
compute_cor <- function(regionlist, wgtlistall, outname = NULL,
                        outputdir = getwd(), ncore = 1) {

  loginfo("Adding R matrix info")
    
  dir.create(file.path(outputdir, paste0(outname, "_LDR")), showWarnings = F)
    
  regionlist_all <- list()
    
  for (b in unique(names(regionlist))){
    regionlist_all[[b]] <- cbind(b, names(regionlist[[b]]))
  }
    
  regionlist_all <- as.data.frame(do.call(rbind, regionlist_all))
  colnames(regionlist_all) <- c("b", "rn")
  #regionlist_all$b <- as.integer(regionlist_all$b)
  #regionlist_all$b <- regionlist_all$b
    
  corelist <- lapply(1:ncore, function(core){njobs <- ceiling(nrow(regionlist_all)/ncore); jobs <- ((core-1)*njobs+1):(core*njobs); jobs[jobs<=nrow(regionlist_all)]})
  names(corelist) <- 1:ncore
    
  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)
    
  outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("tools"), .export=c("read_LD_SNP_file")) %dopar% {
      
    regionlist_core <- regionlist_all[corelist[[core]],]
    outlist_core <- list()  
    b_core <- unique(regionlist_core$b)
      
    for (b in b_core){
      logging::loginfo("Adding R matrix info for chrom %s (core %s)", b, core)  
      regionlist_core_b <- regionlist_core$rn[regionlist_core$b==b]
        
      for (rn in regionlist_core_b){
        outlist_core_region <- list(b=b, rn=rn) 
        R_snp <- lapply(regionlist[[b]][[rn]][["LD_matrix"]], readRDS)
          
        if (length(R_snp)==1){
          R_snp <- unname(R_snp[[1]])
        } else {
          R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
        }
          
        ld_snpinfo <- ctwas:::read_LD_SNP_file(regionlist[[b]][[rn]][["SNP_info"]])
        #sidx <-  match(regionlist[[b]][[rn]][["sid"]], ld_snpinfo$id)
        #outlist_core_region[["sidx"]] <- sidx
        gnames <- regionlist[[b]][[rn]][["gid"]]
        R_snp_gene <- matrix(NA, nrow(R_snp), length(gnames))
        R_gene <- diag(length(gnames))
          
        if (length(gnames) > 0) {
          ldr <- list()
          for (i in 1:length(gnames)){
            gname <- gnames[i]
            wgt <- wgtlistall[[gname]]
            snpnames <- rownames(wgt)
            ld.idx <- match(snpnames, ld_snpinfo$id)
            ldr[[gname]] <- ld.idx
            R.s <- R_snp[ld.idx, ld.idx]
            R_snp_gene[,i] <- sapply(1:nrow(R_snp), function(x){t(wgt)%*%R_snp[ld.idx,x]/sqrt(t(wgt)%*%R.s%*%wgt*R_snp[x,x])})
          }
            
          if (length(gnames) > 1){
            gene_pairs <- combn(length(gnames), 2)
            wgtr <- wgtlistall[gnames]
              
            gene_corrs <- apply(gene_pairs, 2, function(x){t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]/(
              sqrt(t(wgtr[[x[1]]])%*%R_snp[ldr[[x[1]]], ldr[[x[1]]]]%*%wgtr[[x[1]]]) *
                sqrt(t(wgtr[[x[2]]])%*%R_snp[ldr[[x[2]]], ldr[[x[2]]]]%*%wgtr[[x[2]]]))})
            R_gene[t(gene_pairs)] <- gene_corrs
            R_gene[t(gene_pairs[c(2,1),])] <- gene_corrs
          }
        }
    
        #outlist_core_region[["gidx"]] <- regionlist[[b]][[rn]][["gidx"]]
        #outlist_core_region[["gid"]] <- regionlist[[b]][[rn]][["gid"]]
        #outlist_core_region[["LD_matrix"]] <- regionlist[[b]][[rn]][["LD_matrix"]]
        #outlist_core_region[["SNP_info"]] <- regionlist[[b]][[rn]][["SNP_info"]]
        R_sg_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp_gene.RDS"))
        R_g_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_gene.RDS"))
        R_s_file <- file.path(outputdir, paste0(outname, "_LDR"), paste0("chr", b, "_reg", rn, ".R_snp.RDS"))
        R_snp <- R_snp[sidx, sidx, drop = F]
        saveRDS(R_snp, file=R_s_file)
        saveRDS(R_snp_gene, file=R_sg_file)
        saveRDS(R_gene, file=R_g_file)
          
        outlist_core_region[["R_sg_file"]] <- R_sg_file
        outlist_core_region[["R_g_file"]] <- R_g_file
        outlist_core_region[["R_s_file"]] <- R_s_file
        outlist_core[[length(outlist_core)+1]] <- outlist_core_region
      }
    }
    outlist_core
  }

  parallel::stopCluster(cl)
    
  for (i in 1:length(outlist)){
    b <- outlist[[i]][["b"]]
    rn <- outlist[[i]][["rn"]]
    regionlist[[b]][[rn]][["R_sg_file"]] <- outlist[[i]][["R_sg_file"]]
    regionlist[[b]][[rn]][["R_g_file"]] <- outlist[[i]][["R_g_file"]]
    regionlist[[b]][[rn]][["R_s_file"]] <- outlist[[i]][["R_s_file"]]
  }
  return(regionlist)
}
