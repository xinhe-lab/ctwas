
#' Merge region data for cross-boundary genes
#'
#' @param region_data a list of region gene IDs and SNP IDs and associated file names
#'
#' @param boundary_genes a data frame of boundary gene info
#'
#' @param region_info a data frame of region definition and associated LD file names
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param expand TRUE/FALSE. If TRUE, expand merged region_data with full SNPs
#'
#' @param trim_by remove SNPs if the total number of SNPs exceeds limit, options: "random",
#' or "z" (trim SNPs with lower |z|) See parameter `maxSNP` for more information.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param seed seed for random sampling
#'
#' @importFrom logging loginfo
#'
#' @return a list of merged region_data, and a data frame of merged region info
#'
#' @export
#'
merge_region_data <- function(merged_region_list,
                              region_data,
                              region_info,
                              z_snp,
                              z_gene,
                              expand = TRUE,
                              maxSNP = Inf,
                              trim_by = c("random", "z"),
                              ncore = 1,
                              seed = 99) {

  # Merge region data
  loginfo("Merge region_data for %d merged regions", length(merged_region_list))

  merged_region_data <- list()
  for (i in 1:length(merged_region_list)){
    merged_regions <- merged_region_list[[i]]
    new_region_id <- names(merged_region_list)[i]
    region_data_to_merge <- region_data[merged_regions$region_ids]
    minpos <- min(sapply(region_data_to_merge, "[[", "minpos"))
    maxpos <- max(sapply(region_data_to_merge, "[[", "maxpos"))
    gid <- unique(unlist(lapply(region_data_to_merge, "[[", "gid")))
    sid <- unique(unlist(lapply(region_data_to_merge, "[[", "sid")))
    thin <- min(sapply(region_data_to_merge, "[[", "thin"))
    merged_region_data[[new_region_id]] <- list("region_id" = merged_regions$region_ids,
                                                "chrom" = merged_regions$chrom,
                                                "start" = merged_regions$start,
                                                "stop" = merged_regions$stop,
                                                "n_regions" = merged_regions$n_regions,
                                                "minpos" = minpos,
                                                "maxpos" = maxpos,
                                                "gid" = gid,
                                                "sid" = sid,
                                                "thin" = thin)
  }

  if (isTRUE(expand)){
    loginfo("Expand region_data with full SNPs for %d merged regions", length(merged_region_data))
    merged_region_data <- expand_region_data(merged_region_data,
                                             region_info,
                                             z_snp,
                                             z_gene,
                                             trim_by = "z",
                                             maxSNP = maxSNP,
                                             ncore = ncore)
  }else{
    # trim regions with SNPs more than maxSNP
    merged_region_data <- trim_region_data(merged_region_data, z_snp, trim_by = trim_by, maxSNP = maxSNP, seed = seed)

    # add z-scores to region_data
    loginfo("Add z-scores to region_data...")
    merged_region_data <- add_z_to_region_data(merged_region_data, z_snp, z_gene, ncore = ncore)
  }

  return(merged_region_data)

}

#' Identify overlapping regions
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

#' Get a list of regions to be merged
#'
#' @param boundary_genes a data frame of boundary gene info
#'
#' @importFrom logging loginfo
#'
#' @return a list of merged_region_list, and a data frame of merged region info
#'
#' @export
#'
get_merged_regions <- function(boundary_genes){
  # Identify overlapping regions and get a list of regions to be merged
  loginfo("Identify overlapping regions and get a list of regions to be merged.")

  # Identify overlapping regions
  boundary_genes <- label_overlapping_regions(boundary_genes)

  # For each new region, get the region IDs for the regions to be merged
  merged_region_list <- list()
  merge_labels <- unique(boundary_genes$merge_label)
  for( label in merge_labels){
    df <- boundary_genes[boundary_genes$merge_label == label,,drop=F]
    new_region_chrom <- df$chrom[1]
    new_region_start <- min(df$region_start)
    new_region_stop <- max(df$region_stop)
    merged_region_ids <- unique(unlist(strsplit(df[df$merge_label == label, "region_id"], split = ";")))
    new_region_id <- paste0(new_region_chrom, ":", new_region_start, "-", new_region_stop)
    merged_region_list[[new_region_id]] <- list(chrom = new_region_chrom,
                                            start = new_region_start,
                                            stop = new_region_stop,
                                            region_ids = merged_region_ids)
  }

  merged_region_info <- data.frame(chrom = sapply(merged_region_list, "[[", "chrom"),
                                   start = sapply(merged_region_list, "[[", "start"),
                                   stop = sapply(merged_region_list, "[[", "stop"),
                                   region_id = names(merged_region_list))

  for(i in 1:nrow(merged_region_info)){
    merged_region_info[i, "n_merged_regions"] <- length(merged_region_list[[i]]$region_ids)
    merged_region_info[i, "merged_region_ids"] <- paste(merged_region_list[[i]]$region_ids, collapse = ";")
    region_idx <- match(merged_region_list[[i]]$region_ids, region_info$region_id)
    merged_region_info[i, "LD_matrix"] <- paste(region_info[region_idx, "LD_matrix"], collapse = ";")
    merged_region_info[i, "SNP_info"] <- paste(region_info[region_idx, "SNP_info"], collapse = ";")
  }
  rownames(merged_region_info) <- NULL

  return(list(merged_region_list = merged_region_list,
              merged_region_info = merged_region_info))
}

merge_finemap_regions <- function(region_data,
                                  region_info,
                                  boundary_genes,
                                  finemap_res,
                                  z_snp,
                                  z_gene,
                                  weights,
                                  maxSNP = Inf,
                                  expand = TRUE,
                                  group_prior = NULL,
                                  group_prior_var = NULL,
                                  L = 5,
                                  force_compute_cor = FALSE,
                                  save_cor = FALSE,
                                  cor_dir = getwd(),
                                  ncore = 1,
                                  verbose = FALSE){

  high_PIP_genes <- unique(finemap_res[finemap_res$type != "SNP" & finemap_res$susie_pip >= 0.5, "id"])
  high_PIP_boundary_genes <- boundary_genes[boundary_genes$id %in% high_PIP_genes, , drop=FALSE]
  if (nrow(high_PIP_boundary_genes) > 0){
    res <- merge_region_data(region_data,
                             high_PIP_boundary_genes,
                             region_info,
                             z_snp,
                             z_gene,
                             expand = expand,
                             maxSNP = maxSNP)
    merged_region_data <- res$merged_region_data
  }

  finemap_merged_regions_res <- finemap_regions(merged_region_data,
                                                region_info,
                                                weights,
                                                group_prior = group_prior,
                                                group_prior_var = group_prior_var,
                                                L = 5,
                                                force_compute_cor = force_compute_cor,
                                                save_cor = save_cor,
                                                cor_dir = cor_dir,
                                                ncore = ncore,
                                                verbose = verbose)

  return(finemap_merged_regions_res)
}
