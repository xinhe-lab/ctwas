
#' Detect LD mismatches using SuSiE RSS
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#' @param ld_Rinfo a vector of paths to the variant information for all LD matrices
#' @param gwas_n integer, GWAS sample size
#' @param ncore integer, number of cores for parallel computing.
#' @param p_diff_thresh numeric, p-value threshold for identifying problematic SNPs
#' with significant difference between observed z-scores and estimated values
#'
#' @importFrom logging addHandler loginfo
#' @importFrom foreach %dopar% foreach
#'
#' @export
#'
detect_ld_mismatch_susie <- function(z_snp,
                                     region_info,
                                     gwas_n = NULL,
                                     ncore = 1,
                                     p_diff_thresh = 5e-8){

  region_ids <- region_info$region_id
  loginfo("Run LD mismatch diagnosis in %d regions", length(region_ids))

  nregions <- length(region_ids)
  corelist <- lapply(1:ncore, function(core){
    njobs <- ceiling(nregions/ncore);
    jobs <- ((core-1)*njobs+1):(core*njobs);
    jobs[jobs<=nregions]
  })
  names(corelist) <- 1:ncore

  cl <- parallel::makeCluster(ncore, outfile = "")
  doParallel::registerDoParallel(cl)

  outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas", "stats")) %dopar% {

    region_ids_core <- region_ids[corelist[[core]]]

    outlist_core <- list()
    for(region_id in region_ids_core) {

      # Load reference LD matrix and SNP info in the region
      R_snp <- load_LD(region_info$LD_matrix[region_info$region_id == region_id])
      ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info[region_info$region_id == region_id])

      # Match GWAS sumstats with LD reference files. Only keep variants included in LD reference.
      z_snp.region <- z_snp[z_snp$id %in% ld_snpinfo$id,]
      R_snp.idx <- match(z_snp.region$id, ld_snpinfo$id)
      R_snp.region <- R_snp[R_snp.idx, R_snp.idx]
      stopifnot(nrow(z_snp.region) == nrow(R_snp.region))

      # # Estimate lambda (consistency) between the z-scores and LD matrix
      # lambda <- estimate_s_rss(z = z.locus$z, R = R.locus, n = gwas_n)

      # Compute expected z-scores based on conditional distribution of z-scores
      condz_dist <- kriging_rss(z = z_snp.region$z, R = R_snp.region, n = gwas_n)$conditional_dist
      condz_dist <- cbind(z_snp.region[,c("id", "A1", "A2")], condz_dist)

      # compute p-values for the significance of z-score difference between observed and estimated values
      condz_dist$p_diff <- pchisq(condz_dist$z_std_diff^2, df = 1, lower.tail=F)

      outlist_core[[as.character(region_id)]] <- condz_dist
    }
    outlist_core
  }
  parallel::stopCluster(cl)
  stopifnot(length(outlist) == length(region_ids))

  # return problematic variants and flipped variants
  condz_dist <- data.table::rbindlist(outlist, idcol = "region_id")
  problematic_snps <- condz_dist$id[which(condz_dist$p_diff < p_diff_thresh)]
  flipped_snps <- condz_dist$id[which(condz_dist$logLR > 2 & abs(condz_dist$z) > 2)]

  return(list(condz_dist = condz_dist,
              problematic_snps = problematic_snps,
              flipped_snps = flipped_snps))
}


#' Get regions with problematic high PIP SNPs or genes
#'
#' @param problematic_snps a character vector of problematic SNP rsIDs
#' @param highPIP_finemap_res a data frame of the cTWAS finemapping result with high PIP SNPs and genes
#' @param weights weights
#'
#' @return a character vector of region ids with problematic high PIP SNPs or genes
#'
#' @importFrom logging loginfo
#'
#' @export
select_problematic_regions <- function(problematic_snps, highPIP_finemap_res, weights){

  if (length(problematic_snps) == 0) {
    loginfo('No problematic SNPs')
    problematic_ids <- NULL
    problematic_region_ids <- NULL
  }else{
    loginfo('Number of problematic SNPs: %d', length(problematic_snps))

    # find high PIP SNPs that are problematic
    ctwas_highpip_snp_res <- highPIP_finemap_res[highPIP_finemap_res$type == "SNP", ]
    problematic_highpip_snps <- intersect(ctwas_highpip_snp_res$id, problematic_snps)
    loginfo('Number of high PIP SNPs that are problematic: %d', length(problematic_highpip_snps))
    problematic_highpip_snp_region_ids <- unique(highPIP_finemap_res[highPIP_finemap_res$id %in% problematic_highpip_snps,"region_id"])
    loginfo('Number of regions with high PIP SNPs that are problematic: %d', length(problematic_highpip_snp_region_ids))

    # find high PIP genes with problematic SNPs in its weights
    ctwas_highpip_gene_res <- highPIP_finemap_res[highPIP_finemap_res$type != "SNP", ]
    highpip_gids <- ctwas_highpip_gene_res$id
    # subset weights to high PIP genes
    weights_highpip_genes <- weights[highpip_gids]
    # extract snp ids in weights
    problematic_highpip_genes <- c()
    if (length(weights_highpip_genes) > 0){
      for (i in 1:length(weights_highpip_genes)){
        gid <- names(weights_highpip_genes)[i]
        wgt <- weights_highpip_genes[[i]]$wgt
        wgt_snpnames <- rownames(wgt)
        if (any(wgt_snpnames %in% problematic_snps)){
          problematic_highpip_genes <- c(problematic_highpip_genes, gid)
        }
      }
    }
    loginfo('Number of high PIP genes with problematic SNPs in its weights: %d', length(problematic_highpip_genes))
    problematic_highpip_gene_region_ids <- unique(highPIP_finemap_res[highPIP_finemap_res$id %in% problematic_highpip_genes,"region_id"])
    loginfo('Number of regions with high PIP SNPs that are problematic: %d', length(problematic_highpip_gene_region_ids))

    problematic_ids <- unique(c(problematic_highpip_snps, problematic_highpip_genes))

    # get problematic high PIP regions
    problematic_region_ids <- unique(c(problematic_highpip_snp_region_ids, problematic_highpip_gene_region_ids))
    loginfo('Number of problematic regions: %d', length(problematic_region_ids))
  }

  # return problematic region ids
  return(problematic_region_ids = problematic_region_ids)
}

#' Update finemapping result for problematic regions
#'
#' @param old_finemap_res a data frame of old finemapping result
#' @param new_finemap_res a data frame of new finemapping result
#'
#' @export
#'
update_finemap_res <- function(old_finemap_res, new_finemap_res){

  if (!all(colnames(finemap_res) != colnames(new_finemap_res))) {
    stop("columns of finemap_res and new_finemap_res do not match!")
  }

  updated_region_ids <- unique(new_finemap_res$region_id)
  kept_finemap_res <- finemap_res[!finemap_res$region_id %in% updated_region_ids, ]
  updated_finemap_res <- rbind(kept_finemap_res, new_finemap_res)

  return(updated_finemap_res)
}
