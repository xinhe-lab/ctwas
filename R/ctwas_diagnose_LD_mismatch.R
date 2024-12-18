
#' @title Diagnose LD mismatch using SuSiE RSS
#'
#' @param region_ids A vector of region IDs to run diagnosis
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for each of the regions.
#'
#' @param gwas_n integer, GWAS sample size.
#'
#' @param p_diff_thresh numeric, p-value threshold for identifying problematic SNPs
#' with significant difference between observed z-scores and estimated values
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @return a list of problematic SNPs, flipped SNPs,
#' and test statistics from susie's `kriging_rss` function
#'
#' @importFrom logging addHandler loginfo
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#'
#' @export
#'
diagnose_LD_mismatch_susie <- function(region_ids,
                                       z_snp,
                                       LD_map,
                                       gwas_n,
                                       p_diff_thresh = 5e-8,
                                       LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                       LD_loader_fun = NULL,
                                       snpinfo_loader_fun = NULL,
                                       ncore = 1,
                                       logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo("Performing LD mismatch diagnosis for %d regions", length(region_ids))

  if (!all(region_ids %in% LD_map$region_id)){
     stop("region_ids do not match with LD_map$region_id!")
  }

  LD_format <- match.arg(LD_format)

  condz_list <- mclapply_check(region_ids, function(region_id){
    compute_region_condz(region_id, z_snp, LD_map, gwas_n,
                         LD_format = LD_format,
                         LD_loader_fun = LD_loader_fun,
                         snpinfo_loader_fun = snpinfo_loader_fun)
  }, mc.cores = ncore, stop_if_missing = TRUE)

  names(condz_list) <- region_ids
  condz_stats <- rbindlist(condz_list, idcol = "region_id")
  rownames(condz_stats) <- NULL

  # return problematic variants and flipped variants
  problematic_snps <- condz_stats$id[which(condz_stats$p_diff < p_diff_thresh)]
  flipped_snps <- condz_stats$id[which(condz_stats$logLR > 2 & abs(condz_stats$z) > 2)]

  return(list("condz_stats" = condz_stats,
              "problematic_snps" = problematic_snps,
              "flipped_snps" = flipped_snps))
}

# Compute expected z-scores based on conditional distribution of
# z-scores using SuSiE RSS.
#
#' @importFrom stats pchisq
#' @importFrom Matrix bdiag
compute_region_condz <- function(region_id,
                                 z_snp,
                                 LD_map,
                                 gwas_n,
                                 LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                 LD_loader_fun = NULL,
                                 snpinfo_loader_fun = NULL){

  LD_format <- match.arg(LD_format)

  # load LD matrix
  LD_matrix_files <- unlist(strsplit(LD_map$LD_file[LD_map$region_id == region_id], split = ","))
  stopifnot(all(file.exists(LD_matrix_files)))

  if (length(LD_matrix_files) > 1) {
    R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader_fun = LD_loader_fun)
    R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
  } else {
    R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader_fun = LD_loader_fun)
  }

  # load SNP info of the region
  SNP_info_files <- unlist(strsplit(LD_map$SNP_file[LD_map$region_id == region_id], split = ","))
  stopifnot(all(file.exists(SNP_info_files)))
  snpinfo <- read_snp_info_files(SNP_info_files, snpinfo_loader_fun = snpinfo_loader_fun)

  # Match GWAS sumstats with LD reference files. Only keep variants included in LD reference.
  region_z_snp <- z_snp[z_snp$id %in% snpinfo$id,]
  region_z <- region_z_snp$z
  sidx <- match(region_z_snp$id, snpinfo$id)
  region_R <- R_snp[sidx, sidx]

  # Compute expected z-scores based on conditional distribution of z-scores
  condz_stats <- kriging_rss(z = region_z, R = region_R, n = gwas_n)$conditional_dist
  condz_stats <- cbind(region_z_snp[,c("id", "A1", "A2")], condz_stats)

  # compute p-values for the significance of z-score difference between observed and estimated values
  condz_stats$p_diff <- pchisq(condz_stats$z_std_diff^2, df = 1, lower.tail=FALSE)

  return(condz_stats)
}

#' @title Gets problematic genes from problematic SNPs
#'
#' @param problematic_snps a character vector of problematic SNP rsIDs
#'
#' @param weights a list of weights
#'
#' @param finemap_res a data frame of cTWAS finemapping results
#'
#' @param pip_thresh cutoff of gene PIP to select genes
#'
#' @param z_thresh cutoff of abs(z-scores) to select genes
#'
#' @return a vector of problematic genes
#'
#' @importFrom logging loginfo
#'
#' @export
#'
get_problematic_genes <- function(problematic_snps,
                                  weights,
                                  finemap_res,
                                  pip_thresh = 0.5,
                                  z_thresh = NULL){

  if (length(problematic_snps) == 0) {
    loginfo('No problematic SNPs')
    problematic_genes <- NULL
  }else{
    loginfo('Number of problematic SNPs: %d', length(problematic_snps))

    finemap_gene_res <- finemap_res[finemap_res$group!="SNP",]

    # find high PIP genes with problematic SNPs in its weights
    if (!is.null(pip_thresh)) {
      high_pip_gids <- finemap_gene_res$id[finemap_gene_res$susie_pip > pip_thresh]
    } else {
      high_pip_gids <- NULL
    }

    # find high |z| genes with problematic SNPs in its weights
    if (!is.null(z_thresh)) {
      large_z_gids <- finemap_gene_res$id[abs(finemap_gene_res$z) > z_thresh]
    } else {
      large_z_gids <- NULL
    }

    # check both high PIP genes and high |z| genes
    selected_gids <- unique(c(high_pip_gids, large_z_gids))

    selected_weights <- weights[selected_gids]
    # extract SNP ids in weights, and find genes with problematic SNPs in weights
    problematic_genes <- NULL
    if (length(selected_weights) > 0){
      for (i in 1:length(selected_weights)){
        gid <- names(selected_weights)[i]
        wgt <- selected_weights[[i]]$wgt
        wgt_snp_ids <- rownames(wgt)
        if (any(wgt_snp_ids %in% problematic_snps)){
          problematic_genes <- c(problematic_genes, gid)
        }
      }
    }
    loginfo('Number of problematic genes: %d', length(problematic_genes))
  }

  return(problematic_genes)
}


#' Updates cTWAS finemapping result for selected regions
#'
#' @param finemap_res a data frame of original finemapping result.
#' @param susie_alpha_res a data frame of original susie alpha result.
#' @param new_finemap_res a data frame of new finemapping result.
#' @param new_susie_alpha_res a data frame of new susie alpha result.
#' @param updated_region_ids a vector of region ids to be updated.
#'
#' @return a list with updated cTWAS finemapping result.
#'
#' @export
update_finemap_res <- function(finemap_res,
                               susie_alpha_res,
                               new_finemap_res,
                               new_susie_alpha_res,
                               updated_region_ids){

  if (!all(colnames(finemap_res) == colnames(new_finemap_res)))
    stop("columns of finemap_res and new_finemap_res do not match!")

  if (!all(colnames(susie_alpha_res) == colnames(new_susie_alpha_res)))
    stop("columns of susie_alpha_res and new_susie_alpha_res do not match!")

  if (missing(updated_region_ids)){
    updated_region_ids <- unique(new_finemap_res$region_id)
  }

  kept_finemap_res <- finemap_res[!finemap_res$region_id %in% updated_region_ids, ]
  new_finemap_res <- new_finemap_res[new_finemap_res$region_id %in% updated_region_ids, ]
  finemap_res <- rbind(kept_finemap_res, new_finemap_res)

  kept_susie_alpha_res <- susie_alpha_res[!susie_alpha_res$region_id %in% updated_region_ids, ]
  new_susie_alpha_res <- new_susie_alpha_res[new_susie_alpha_res$region_id %in% updated_region_ids, ]
  susie_alpha_res <- rbind(kept_susie_alpha_res, new_susie_alpha_res)

  return(list("finemap_res" = finemap_res,
              "susie_alpha_res" = susie_alpha_res))
}


