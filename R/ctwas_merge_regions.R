#' @title Merges regions with boundary genes and finemaps merged regions
#'
#' @param boundary_genes a data frame of boundary genes
#'
#' @param region_data a list of original region_data
#'
#' @param region_info a data frame of region definitions
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param weights a list of weights for each gene
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions.
#'
#' @param snp_info a list of SNP info data frames for LD reference.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param expand TRUE/FALSE. If TRUE, expand merged region_data with full SNPs
#'
#' @param L the number of effects for susie during the fine mapping steps
#'
#' @param group_prior a vector of two prior inclusion probabilities for SNPs and genes.
#'
#' @param group_prior_var a vector of two prior variances for SNPs and gene effects.
#'
#' @param use_null_weight TRUE/FALSE. If TRUE, allow for a probability of no effect in susie
#'
#' @param coverage A number between 0 and 1 specifying the \dQuote{coverage} of the estimated confidence sets
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param force_compute_cor TRUE/FALSE. If TRUE, force computing correlation (R) matrices
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices
#'
#' @param cor_dir a string, the directory to store correlation (R) matrices
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param ncore The number of cores used to parallelize computation over regions
#'
#' @param verbose TRUE/FALSE. If TRUE, print detail messages
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a list of merged region data, merged region info,
#'   merged region IDs,and finemapping results for merged regions.
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @export
#'
merge_finemap_regions <- function(boundary_genes,
                                  region_data,
                                  region_info,
                                  z_snp,
                                  z_gene,
                                  weights,
                                  use_LD = TRUE,
                                  LD_info = NULL,
                                  snp_info = NULL,
                                  maxSNP = Inf,
                                  expand = TRUE,
                                  L = 5,
                                  group_prior = NULL,
                                  group_prior_var = NULL,
                                  use_null_weight = TRUE,
                                  coverage = 0.95,
                                  min_abs_corr = 0.5,
                                  force_compute_cor = FALSE,
                                  save_cor = FALSE,
                                  cor_dir = getwd(),
                                  LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                                  LD_loader_fun,
                                  ncore = 1,
                                  verbose = FALSE,
                                  ...){

  if (nrow(boundary_genes) == 0) {
    loginfo("No regions to merge")
  } else {
    loginfo("%d boundary genes to merge", nrow(boundary_genes))
    res <- merge_region_data(boundary_genes,
                             region_data,
                             region_info = region_info,
                             use_LD = use_LD,
                             LD_info = LD_info,
                             snp_info = snp_info,
                             z_snp = z_snp,
                             z_gene = z_gene,
                             expand = TRUE,
                             maxSNP = maxSNP)
    merged_region_data <- res$merged_region_data
    merged_region_info <- res$merged_region_info
    merged_LD_info <- res$merged_LD_info
    merged_snp_info <- res$merged_snp_info
    merge_region_id_list <- res$merge_region_id_list

    finemap_merged_regions_res <- finemap_regions(merged_region_data,
                                                  use_LD = use_LD,
                                                  LD_info = merged_LD_info,
                                                  snp_info = merged_snp_info,
                                                  weights = weights,
                                                  group_prior = group_prior,
                                                  group_prior_var = group_prior_var,
                                                  L = L,
                                                  force_compute_cor = force_compute_cor,
                                                  save_cor = save_cor,
                                                  cor_dir = cor_dir,
                                                  LD_format = LD_format,
                                                  LD_loader_fun = LD_loader_fun,
                                                  ncore = ncore,
                                                  verbose = verbose,
                                                  ...)

    return(list(finemap_merged_regions_res = finemap_merged_regions_res,
                merged_region_data = merged_region_data,
                merged_region_info = merged_region_info,
                merged_LD_info = merged_LD_info,
                merged_snp_info = merged_snp_info,
                merge_region_id_list = merge_region_id_list))
  }
}


# Identify overlapping regions
label_overlapping_regions <- function(boundary_genes) {

  boundary_genes <- boundary_genes[,c("id", "chrom", "region_start", "region_stop", "region_id", "n_regions")]
  boundary_genes <- boundary_genes[boundary_genes$n_regions > 1, ]

  boundary_genes_list <- list()
  for(b in unique(boundary_genes$chrom)){
    df <- boundary_genes[boundary_genes$chrom==b, ]

    # Sort the dataframe by 'start' and then by 'stop'
    df <- df[with(df, order(region_start, region_stop)),]

    # Initialize the label counter and assign the first label
    merge_idx <- 1
    df$merge_label <- paste0(b, ":", merge_idx)

    if (nrow(df) > 1) {
      # Loop through the rows starting from the second row
      for (i in 2:nrow(df)) {
        # Check if the current row overlaps with the previous row
        if (df$region_start[i] <= df$region_stop[i-1]) {
          # merge with previous region
          df$merge_label[i] <- paste0(b, ":", merge_idx)
        } else {
          # do not merge, create a new merge_label
          merge_idx <- merge_idx + 1
          df$merge_label[i] <- paste0(b, ":", merge_idx)
        }
      }
    }

    boundary_genes_list[[as.character(b)]] <- df
  }

  boundary_genes <- do.call(rbind, boundary_genes_list)
  rownames(boundary_genes) <- NULL
  return(boundary_genes)
}

# Get merged region info
#' @importFrom logging loginfo
get_merged_region_info <- function(boundary_genes, region_info, use_LD = TRUE, LD_info, snp_info){

  # Identify overlapping regions
  boundary_genes <- label_overlapping_regions(boundary_genes)

  # For each new region, get the region IDs for the regions to be merged
  merge_labels <- unique(boundary_genes$merge_label)

  merged_region_info <- data.frame()
  merged_snp_info <- list()
  merge_region_id_list <- list()
  merged_LD_info <- NULL
  for(merge_label in merge_labels){
    df <- boundary_genes[boundary_genes$merge_label == merge_label,,drop=F]
    new_region_chrom <- df$chrom[1]
    new_region_start <- min(df$region_start)
    new_region_stop <- max(df$region_stop)
    new_region_id <- paste0(new_region_chrom, "_", new_region_start, "_", new_region_stop)

    merged_region_ids <- unique(unlist(strsplit(df[df$merge_label == merge_label, "region_id"], split = ";")))

    merged_region_info <- rbind(merged_region_info,
                                data.frame(chrom = new_region_chrom,
                                           start = new_region_start,
                                           stop = new_region_stop,
                                           region_id = new_region_id))

    merged_snp_info[[new_region_id]] <-  do.call(rbind, snp_info[merged_region_ids])

    merge_region_id_list[[new_region_id]] <- list(region_id = new_region_id,
                                                  merged_region_ids = merged_region_ids)

    if (use_LD) {
      merged_region_idx <- match(merged_region_ids, LD_info$region_id)
      merged_LD_matrix_files <- paste(LD_info$LD_matrix[merged_region_idx], collapse = ";")
      merged_LD_info <- rbind(merged_LD_info,
                              data.frame(region_id = new_region_id, LD_matrix = merged_LD_matrix_files))
    }
  }

  loginfo("Merged %d boundary genes into %d regions", nrow(boundary_genes), nrow(merged_region_info))

  return(list(merged_region_info = merged_region_info,
              merged_LD_info = merged_LD_info,
              merged_snp_info = merged_snp_info,
              merge_region_id_list = merge_region_id_list))
}

#' @title Merges region data for cross-boundary genes
#'
#' @param boundary_genes a data frame of boundary gene info
#'
#' @param region_data a list of original region_data
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions.
#'
#' @param snp_info a list of SNP info data frames for LD reference.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param expand TRUE/FALSE. If TRUE, expand merged region_data with full SNPs
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit, options: "random",
#' or "z" (trim SNPs with lower |z|) See parameter `maxSNP` for more information.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @return a list of merged region data, merged region info, LD info, snp info,
#' and merged region IDs.
#'
#' @importFrom logging loginfo
#'
#' @export
#'
merge_region_data <- function(boundary_genes,
                              region_data,
                              region_info,
                              use_LD = TRUE,
                              LD_info,
                              snp_info,
                              z_snp,
                              z_gene,
                              expand = TRUE,
                              maxSNP = Inf,
                              trim_by = c("random", "z"),
                              ncore = 1,
                              seed = 99) {

  # Identify overlapping regions and get a list of regions to be merged
  loginfo("Identify overlapping regions and get region_info for merged regions.")
  res <- get_merged_region_info(boundary_genes, region_info, use_LD, LD_info, snp_info)
  merged_region_info <- res$merged_region_info
  merged_LD_info <- res$merged_LD_info
  merged_snp_info <- res$merged_snp_info
  merge_region_id_list <- res$merge_region_id_list

  # Merge region data
  loginfo("Merging region_data ...")
  merged_region_data <- list()
  for (i in 1:nrow(merged_region_info)){
    merged_regioninfo <- merged_region_info[i,]
    region_id <- merged_regioninfo$region_id
    chrom <- merged_regioninfo$chrom
    start <- merged_regioninfo$start
    stop <- merged_regioninfo$stop
    merged_region_ids <- merge_region_id_list[[region_id]][["merged_region_ids"]]
    # merge gids and sids from the old region_data
    gid <- unique(unlist(lapply(region_data[merged_region_ids], "[[", "gid")))
    sid <- unique(unlist(lapply(region_data[merged_region_ids], "[[", "sid")))
    minpos <- min(sapply(region_data[merged_region_ids], "[[", "minpos"))
    maxpos <- max(sapply(region_data[merged_region_ids], "[[", "maxpos"))
    thin <- min(sapply(region_data[merged_region_ids], "[[", "thin"))
    merged_region_data[[region_id]] <- list("region_id" = region_id,
                                            "chrom" = chrom,
                                            "start" = start,
                                            "stop" = stop,
                                            "minpos" = minpos,
                                            "maxpos" = maxpos,
                                            "gid" = gid,
                                            "sid" = sid,
                                            "thin" = thin)
  }
  loginfo("%d regions in merged_region_data", length(merged_region_data))

  if (expand) {
    # expand with full SNPs
    merged_region_data <- expand_region_data(merged_region_data,
                                             merged_snp_info,
                                             z_snp,
                                             z_gene,
                                             trim_by = "z",
                                             maxSNP = maxSNP,
                                             ncore = ncore)
  }

  # trim regions with SNPs more than maxSNP
  merged_region_data <- trim_region_data(merged_region_data, z_snp,
                                         trim_by = trim_by,
                                         maxSNP = maxSNP,
                                         seed = seed)

  # add z-scores to region_data
  merged_region_data <- add_z_to_region_data(merged_region_data,
                                             z_snp, z_gene, ncore = ncore)

  return(list(merged_region_data = merged_region_data,
              merged_region_info = merged_region_info,
              merged_LD_info = merged_LD_info,
              merged_snp_info = merged_snp_info,
              merge_region_id_list = merge_region_id_list))

}

#' @title Updates finemapping result for merged regions
#'
#' @param finemap_res a data frame of finemapping result
#' @param finemap_merged_regions_res a data frame of finemapping result for merged regions
#' @param merge_region_id_list a list with new region IDs and merged region IDs
#'
#' @export
#'
update_merged_regions_finemap_res <- function(finemap_res,
                                              finemap_merged_regions_res,
                                              merge_region_id_list){
  for(region_id in names(merge_region_id_list)){
    merged_region_ids <- merge_region_id_list[[region_id]][["merged_region_ids"]]
    finemap_res[finemap_res$region_id %in% merged_region_ids, ] <-
      finemap_merged_regions_res[finemap_merged_regions_res$region_id==region_id, ]
  }
  return(finemap_res)
}
