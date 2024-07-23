
#' @title Diagnose LD mismatches using SuSiE RSS
#'
#' @param z_snp A data frame with two columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param region_ids A vector of region IDs to run diagnosis
#'
#' @param LD_map a data frame with filenames of LD matrices for each of the regions.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param gwas_n integer, GWAS sample size.
#'
#' @param ncore integer, number of cores for parallel computing.
#'
#' @param p_diff_thresh numeric, p-value threshold for identifying problematic SNPs
#' with significant difference between observed z-scores and estimated values
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
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
diagnose_ld_mismatch_susie <- function(z_snp,
                                       region_ids,
                                       LD_map,
                                       snp_map,
                                       gwas_n = NULL,
                                       ncore = 1,
                                       p_diff_thresh = 5e-8,
                                       LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                       LD_loader_fun,
                                       logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  loginfo("Perform LD mismatch diagnosis for %d regions", length(region_ids))
  LD_format <- match.arg(LD_format)

  condz_list <- mclapply_check(region_ids, function(region_id){
    compute_region_condz(region_id, LD_map, snp_map, z_snp, gwas_n,
                         LD_format = LD_format,
                         LD_loader_fun = LD_loader_fun)
  }, mc.cores = ncore)

  names(condz_list) <- region_ids
  condz_stats <- rbindlist(condz_list, idcol = "region_id")
  rownames(condz_stats) <- NULL

  # return problematic variants and flipped variants
  problematic_snps <- condz_stats$id[which(condz_stats$p_diff < p_diff_thresh)]
  flipped_snps <- condz_stats$id[which(condz_stats$logLR > 2 & abs(condz_stats$z) > 2)]

  return(list(condz_stats = condz_stats,
              problematic_snps = problematic_snps,
              flipped_snps = flipped_snps))
}

# Compute expected z-scores based on conditional distribution of
# z-scores using SuSiE RSS.
#
#' @importFrom stats pchisq
#' @importFrom Matrix bdiag
compute_region_condz <- function(region_id, LD_map, snp_map, z_snp, gwas_n,
                                 LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                 LD_loader_fun){

  LD_format <- match.arg(LD_format)

  # load LD matrix
  LD_matrix_files <- unlist(strsplit(LD_map[LD_map$region_id == region_id, "LD_file"], split = ";"))
  stopifnot(all(file.exists(LD_matrix_files)))
  if (length(LD_matrix_files) > 1) {
    R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader_fun = LD_loader_fun)
    R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
  } else {
    R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader_fun = LD_loader_fun)
  }

  # load SNP info
  snpinfo <- do.call(rbind, snp_map[region_id])

  # Match GWAS sumstats with LD reference files. Only keep variants included in LD reference.
  region_z_snp <- z_snp[z_snp$id %in% snpinfo$id,]
  sidx <- match(region_z_snp$id, snpinfo$id)
  region_R_snp <- R_snp[sidx, sidx]

  # # Estimate lambda (consistency) between the z-scores and LD matrix
  # lambda <- estimate_s_rss(z = z.locus$z, R = R.locus, n = gwas_n)

  # Compute expected z-scores based on conditional distribution of z-scores
  condz_stats <- kriging_rss(z = region_z_snp$z, R = region_R_snp, n = gwas_n)$conditional_dist
  condz_stats <- cbind(region_z_snp[,c("id", "A1", "A2")], condz_stats)

  # compute p-values for the significance of z-score difference between observed and estimated values
  condz_stats$p_diff <- pchisq(condz_stats$z_std_diff^2, df = 1, lower.tail=F)

  return(condz_stats)
}

#' @title Gets problematic genes from problematic SNPs
#'
#' @param problematic_snps a character vector of problematic SNP rsIDs.
#'
#' @param weights a list of weights
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
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

# Updates finemapping result
update_finemap_res <- function(finemap_res, new_finemap_res){

  if (!all(colnames(finemap_res) == colnames(new_finemap_res))) {
    stop("columns of finemap_res and new_finemap_res do not match!")
  }

  updated_region_ids <- unique(new_finemap_res$region_id)
  kept_finemap_res <- finemap_res[!finemap_res$region_id %in% updated_region_ids, ]
  updated_finemap_res <- rbind(kept_finemap_res, new_finemap_res)

  return(updated_finemap_res)
}
