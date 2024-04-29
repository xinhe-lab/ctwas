
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
#' @return a list of problematic SNPs, flipped SNPs,
#' and test statistics from susie's `kriging_rss` function
#'
#' @export
#'
detect_ld_mismatch_susie <- function(z_snp,
                                     region_info,
                                     gwas_n = NULL,
                                     ncore = 1,
                                     p_diff_thresh = 5e-8){

  region_ids <- region_info$region_id
  loginfo("Perform LD mismatch diagnosis for %d regions", length(region_ids))

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

      regioninfo <- region_info[region_info$region_id %in% region_id, ]

      # load LD matrix
      LD_matrix_files <- unlist(strsplit(regioninfo$LD_matrix, split = ";"))
      stopifnot(all(file.exists(LD_matrix_files)))
      if (length(LD_matrix_files) == 1){
        R_snp <- load_LD(LD_matrix_files)
      } else {
        R_snp <- lapply(LD_matrix_files, load_LD)
        R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
      }

      # load SNP info
      SNP_info_files <- unlist(strsplit(regioninfo$SNP_info, split = ";"))
      stopifnot(all(file.exists(SNP_info_files)))
      ld_snpinfo <- read_LD_SNP_files(SNP_info_files)

      # Match GWAS sumstats with LD reference files. Only keep variants included in LD reference.
      z_snp.region <- z_snp[z_snp$id %in% ld_snpinfo$id,]
      R_snp.idx <- match(z_snp.region$id, ld_snpinfo$id)
      R_snp.region <- R_snp[R_snp.idx, R_snp.idx]
      stopifnot(nrow(z_snp.region) == nrow(R_snp.region))

      # # Estimate lambda (consistency) between the z-scores and LD matrix
      # lambda <- estimate_s_rss(z = z.locus$z, R = R.locus, n = gwas_n)

      # Compute expected z-scores based on conditional distribution of z-scores
      condz_stats <- kriging_rss(z = z_snp.region$z, R = R_snp.region, n = gwas_n)$conditional_dist
      condz_stats <- cbind(z_snp.region[,c("id", "A1", "A2")], condz_stats)

      # compute p-values for the significance of z-score difference between observed and estimated values
      condz_stats$p_diff <- pchisq(condz_stats$z_std_diff^2, df = 1, lower.tail=F)

      outlist_core[[as.character(region_id)]] <- condz_stats
    }
    outlist_core
  }
  parallel::stopCluster(cl)
  stopifnot(length(outlist) == length(region_ids))

  # return problematic variants and flipped variants
  condz_stats <- data.table::rbindlist(outlist, idcol = "region_id")
  problematic_snps <- condz_stats$id[which(condz_stats$p_diff < p_diff_thresh)]
  flipped_snps <- condz_stats$id[which(condz_stats$logLR > 2 & abs(condz_stats$z) > 2)]

  return(list(condz_stats = condz_stats,
              problematic_snps = problematic_snps,
              flipped_snps = flipped_snps))
}


#' Get problematic genes
#'
#' @param problematic_snps a character vector of problematic SNP rsIDs
#' @param weights a list of weights
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#' @param z_thresh cutoff of abs(z-scores) to select genes with large effect sizes
#'
#' @return a vector of problematic genes
#'
#' @importFrom logging loginfo
#'
#' @export
#'
get_problematic_genes <- function(problematic_snps, weights, z_gene, z_thresh = 3){

  if (length(problematic_snps) == 0) {
    loginfo('No problematic SNPs')
    problematic_genes <- NULL
  }else{
    loginfo('Number of problematic SNPs: %d', length(problematic_snps))

    # find high |z| genes with problematic SNPs in its weights
    selected_gids <- z_gene[abs(z_gene$z) > z_thresh, "id"]
    # subset weights to high |z| genes
    selected_weights <- weights[selected_gids]
    # extract snp ids in weights
    problematic_genes <- c()
    if (length(selected_weights) > 0){
      for (i in 1:length(selected_weights)){
        gid <- names(selected_weights)[i]
        wgt <- selected_weights[[i]]$wgt
        wgt_snpnames <- rownames(wgt)
        if (any(wgt_snpnames %in% problematic_snps)){
          problematic_genes <- c(problematic_genes, gid)
        }
      }
    }
    loginfo('Number of large effect genes with problematic SNPs in weights: %d', length(problematic_genes))
  }

  # return problematic gene ids
  return(problematic_genes)
}

#' Update finemapping result
#'
#' @param finemap_res a data frame of original finemapping result
#' @param new_finemap_res a data frame of new finemapping result
#'
#' @export
#'
update_finemap_res <- function(finemap_res, new_finemap_res){

  if (!all(colnames(finemap_res) == colnames(new_finemap_res))) {
    stop("columns of finemap_res and new_finemap_res do not match!")
  }

  updated_region_ids <- unique(new_finemap_res$region_id)
  kept_finemap_res <- finemap_res[!finemap_res$region_id %in% updated_region_ids, ]
  updated_finemap_res <- rbind(kept_finemap_res, new_finemap_res)

  return(updated_finemap_res)
}
