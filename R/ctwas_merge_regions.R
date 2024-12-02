
#' @title Merges region data for cross-boundary genes
#'
#' @param boundary_genes a data frame of boundary gene info
#'
#' @param region_data a list of original region_data
#'
#' @param region_info a data frame of region definitions
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for the regions.
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param weights a list of preprocessed weights.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param estimate_L If TRUE, estimate L for merged regions.
#'
#' @param expand If TRUE, expand merged region_data with full SNPs
#'
#' @param L the number of effects for susie.
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param verbose If TRUE, print detail messages
#'
#' @param logfile The log filename. If NULL, will print log info on screen.
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a list of merged region data, merged region info, LD_map, snp_map,
#' and merged region IDs.
#'
#' @importFrom logging loginfo
#'
#' @export
#'
merge_region_data <- function(boundary_genes,
                              region_data,
                              region_info,
                              LD_map,
                              snp_map,
                              weights,
                              z_snp,
                              z_gene,
                              estimate_L = TRUE,
                              expand = TRUE,
                              L = 5,
                              maxSNP = Inf,
                              LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
                              LD_loader_fun = NULL,
                              snpinfo_loader_fun = NULL,
                              ncore = 1,
                              verbose = FALSE,
                              logfile = NULL,
                              ...) {

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  # Identify overlapping regions and get a list of regions to be merged
  loginfo("Identify overlapping regions and create merged snp_map and LD_map.")
  res <- create_merged_snp_LD_map(boundary_genes,
                                  region_info = region_info,
                                  snp_map = snp_map,
                                  LD_map = LD_map)
  merged_region_info <- res$merged_region_info
  merged_LD_map <- res$merged_LD_map
  merged_snp_map <- res$merged_snp_map
  merged_region_id_map <- res$merged_region_id_map
  rm(res)

  # Merge region data
  loginfo("Merging region data ...")
  merged_region_data <- list()
  for (i in 1:nrow(merged_region_info)){
    merged_regioninfo <- merged_region_info[i,]
    new_region_id <- merged_regioninfo$region_id

    old_region_ids <- merged_region_id_map$old_region_ids[merged_region_id_map$region_id == new_region_id]
    old_region_ids <- unlist(strsplit(old_region_ids, ","))

    new_chrom <- merged_regioninfo$chrom
    new_start <- merged_regioninfo$start
    new_stop <- merged_regioninfo$stop

    # merge gids and sids from the old region_data
    new_gid <- as.character(unlist(lapply(region_data[old_region_ids], "[[", "gid")))
    new_sid <- as.character(unlist(lapply(region_data[old_region_ids], "[[", "sid")))
    new_minpos <- min(sapply(region_data[old_region_ids], "[[", "minpos"))
    new_maxpos <- max(sapply(region_data[old_region_ids], "[[", "maxpos"))
    new_thin <- min(sapply(region_data[old_region_ids], "[[", "thin"))
    new_groups <- unique(as.character(sapply(region_data[old_region_ids], "[[", "groups")))
    merged_region_data[[new_region_id]] <- list("region_id" = new_region_id,
                                                "chrom" = new_chrom,
                                                "start" = new_start,
                                                "stop" = new_stop,
                                                "minpos" = new_minpos,
                                                "maxpos" = new_maxpos,
                                                "thin" = new_thin,
                                                "gid" = new_gid,
                                                "sid" = new_sid)
  }

  merged_region_data <- update_region_z(merged_region_data, z_snp, z_gene, ncore = ncore)

  loginfo("%d regions in merged_region_data", length(merged_region_data))

  if (estimate_L) {
    loginfo("Estimating L ...")
    merged_region_L <- estimate_region_L(region_data = merged_region_data,
                                         LD_map = merged_LD_map,
                                         weights = weights,
                                         init_L = L,
                                         LD_format = LD_format,
                                         LD_loader_fun = LD_loader_fun,
                                         snpinfo_loader_fun = snpinfo_loader_fun,
                                         ncore = ncore,
                                         verbose = verbose,
                                         ...)
    merged_region_L[merged_region_L == 0] <- 1
  } else {
    merged_region_L <- L
  }

  if (expand) {
    merged_region_data <- expand_region_data(merged_region_data,
                                             merged_snp_map,
                                             z_snp = z_snp,
                                             maxSNP = maxSNP,
                                             ncore = ncore)
  }

  return(list("merged_region_data" = merged_region_data,
              "merged_region_info" = merged_region_info,
              "merged_LD_map" = merged_LD_map,
              "merged_snp_map" = merged_snp_map,
              "merged_region_id_map" = merged_region_id_map,
              "merged_region_L" = merged_region_L))
}

#' @title Merges region data for cross-boundary genes without using LD
#'
#' @param boundary_genes a data frame of boundary gene info
#'
#' @param region_data a list of original region_data
#'
#' @param region_info a data frame of region definitions
#'
#' @param snp_map a list of data frames with SNP-to-region map for the reference.
#'
#' @param z_snp A data frame with columns: "id", "z", giving the z-scores for SNPs.
#'
#' @param z_gene A data frame with columns: "id", "z", giving the z-scores for genes.
#'
#' @param expand If TRUE, expand merged region_data with full SNPs
#'
#' @param maxSNP Inf or integer. Maximum number of SNPs in a region. Default is
#' Inf, no limit. This can be useful if there are many SNPs in a region and you don't
#' have enough memory to run the program.
#'
#' @param ncore The number of cores used to parallelize susie over regions
#'
#' @param verbose If TRUE, print detail messages
#'
#' @param logfile The log filename. If NULL, will print log info on screen.
#'
#' @param ... Additional arguments of \code{susie_rss}.
#'
#' @return a list of merged region data, merged region info, snp_map,
#' and merged region IDs.
#'
#' @importFrom logging loginfo
#'
#' @export
#'
merge_region_data_noLD <- function(boundary_genes,
                                   region_data,
                                   region_info,
                                   snp_map,
                                   z_snp,
                                   z_gene,
                                   expand = TRUE,
                                   maxSNP = Inf,
                                   ncore = 1,
                                   verbose = FALSE,
                                   logfile = NULL,
                                   ...) {

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  # Identify overlapping regions and get a list of regions to be merged
  loginfo("Identify overlapping regions and create merged snp_map.")
  res <- create_merged_snp_map(boundary_genes, region_info, snp_map)
  merged_region_info <- res$merged_region_info
  merged_snp_map <- res$merged_snp_map
  merged_region_id_map <- res$merged_region_id_map
  rm(res)

  # Merge region data
  loginfo("Merging region data ...")
  merged_region_data <- list()
  for (i in 1:nrow(merged_region_info)){
    merged_regioninfo <- merged_region_info[i,]
    new_region_id <- merged_regioninfo$region_id

    old_region_ids <- merged_region_id_map$old_region_ids[merged_region_id_map$region_id == new_region_id]
    old_region_ids <- unlist(strsplit(old_region_ids, ","))

    new_chrom <- merged_regioninfo$chrom
    new_start <- merged_regioninfo$start
    new_stop <- merged_regioninfo$stop

    # merge gids and sids from the old region_data
    new_gid <- as.character(unlist(lapply(region_data[old_region_ids], "[[", "gid")))
    new_sid <- as.character(unlist(lapply(region_data[old_region_ids], "[[", "sid")))
    new_minpos <- min(sapply(region_data[old_region_ids], "[[", "minpos"))
    new_maxpos <- max(sapply(region_data[old_region_ids], "[[", "maxpos"))
    new_thin <- min(sapply(region_data[old_region_ids], "[[", "thin"))
    new_groups <- unique(as.character(sapply(region_data[old_region_ids], "[[", "groups")))
    merged_region_data[[new_region_id]] <- list("region_id" = new_region_id,
                                                "chrom" = new_chrom,
                                                "start" = new_start,
                                                "stop" = new_stop,
                                                "minpos" = new_minpos,
                                                "maxpos" = new_maxpos,
                                                "thin" = new_thin,
                                                "gid" = new_gid,
                                                "sid" = new_sid)
  }

  merged_region_data <- update_region_z(merged_region_data, z_snp, z_gene, ncore = ncore)

  loginfo("%d regions in merged_region_data", length(merged_region_data))

  if (expand) {
    merged_region_data <- expand_region_data(merged_region_data,
                                             merged_snp_map,
                                             z_snp = z_snp,
                                             maxSNP = maxSNP,
                                             ncore = ncore)
  }

  return(list("merged_region_data" = merged_region_data,
              "merged_region_info" = merged_region_info,
              "merged_snp_map" = merged_snp_map,
              "merged_region_id_map" = merged_region_id_map))
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

# Get merged region_info, snp_map and LD_map
#' @importFrom logging loginfo
create_merged_snp_LD_map <- function(boundary_genes,
                                     region_info,
                                     snp_map,
                                     LD_map){
  # Identify overlapping regions
  boundary_genes <- label_overlapping_regions(boundary_genes)

  # For each new region, get the region IDs for the regions to be merged
  merge_labels <- unique(boundary_genes$merge_label)

  merged_region_info <- data.frame()
  merged_snp_map <- list()
  merged_LD_map <- data.frame()
  merged_region_id_map <- data.frame()

  for(merge_label in merge_labels){
    df <- boundary_genes[boundary_genes$merge_label == merge_label,,drop=F]
    new_region_chrom <- df$chrom[1]
    new_region_start <- min(df$region_start)
    new_region_stop <- max(df$region_stop)
    new_region_id <- paste0(new_region_chrom, "_", new_region_start, "_", new_region_stop)

    old_region_ids <- unique(unlist(strsplit(df[df$merge_label == merge_label, "region_id"], split = ",")))

    merged_region_info <- rbind(merged_region_info,
                                data.frame(chrom = new_region_chrom,
                                           start = new_region_start,
                                           stop = new_region_stop,
                                           region_id = new_region_id))

    old_LD_files <- paste(LD_map$LD_file[match(old_region_ids, LD_map$region_id)], collapse = ",")
    old_SNP_files <- paste(LD_map$SNP_file[match(old_region_ids, LD_map$region_id)], collapse = ",")
    merged_LD_map <- rbind(merged_LD_map,
                           data.frame(region_id = new_region_id,
                                      LD_file = old_LD_files,
                                      SNP_file = old_SNP_files))

    merged_snp_map[[new_region_id]] <-  do.call(rbind, snp_map[old_region_ids])

    merged_region_id_map <- rbind(merged_region_id_map,
                                  data.frame(region_id = new_region_id,
                                             old_region_ids = paste(old_region_ids, collapse = ",")))

  }

  loginfo("Merge %d boundary genes into %d regions", nrow(boundary_genes), nrow(merged_region_info))

  return(list("merged_region_info" = merged_region_info,
              "merged_LD_map" = merged_LD_map,
              "merged_snp_map" = merged_snp_map,
              "merged_region_id_map" = merged_region_id_map))
}

# Get merged region_info, and snp_map
#' @importFrom logging loginfo
create_merged_snp_map <- function(boundary_genes,
                                  region_info,
                                  snp_map){
  # Identify overlapping regions
  boundary_genes <- label_overlapping_regions(boundary_genes)

  # For each new region, get the region IDs for the regions to be merged
  merge_labels <- unique(boundary_genes$merge_label)

  merged_region_info <- data.frame()
  merged_snp_map <- list()
  merged_region_id_map <- data.frame()
  for(merge_label in merge_labels){
    df <- boundary_genes[boundary_genes$merge_label == merge_label,,drop=F]
    new_region_chrom <- df$chrom[1]
    new_region_start <- min(df$region_start)
    new_region_stop <- max(df$region_stop)
    new_region_id <- paste0(new_region_chrom, "_", new_region_start, "_", new_region_stop)

    old_region_ids <- unique(unlist(strsplit(df[df$merge_label == merge_label, "region_id"], split = ",")))

    merged_region_info <- rbind(merged_region_info,
                                data.frame(chrom = new_region_chrom,
                                           start = new_region_start,
                                           stop = new_region_stop,
                                           region_id = new_region_id))

    merged_snp_map[[new_region_id]] <-  do.call(rbind, snp_map[old_region_ids])

    merged_region_id_map <- rbind(merged_region_id_map,
                                  data.frame(region_id = new_region_id,
                                             old_region_ids = paste(old_region_ids, collapse = ",")))

  }

  loginfo("Merge %d boundary genes into %d regions", nrow(boundary_genes), nrow(merged_region_info))

  return(list("merged_region_info" = merged_region_info,
              "merged_snp_map" = merged_snp_map,
              "merged_region_id_map" = merged_region_id_map))
}



#' Updates cTWAS finemapping result for merged regions
#'
#' @param finemap_res a data frame of original finemapping result.
#' @param susie_alpha_res a data frame of original susie alpha result.
#' @param merged_region_finemap_res a data frame of finemapping result for merged regions.
#' @param merged_region_susie_alpha_res a data frame of susie alpha result for merged regions.
#' @param merged_region_id_map a data frame of new region IDs and original regions IDs.
#'
#' @return a list with updated cTWAS finemapping result.
#'
#' @export
#'
update_merged_region_finemap_res <- function(finemap_res,
                                             susie_alpha_res,
                                             merged_region_finemap_res,
                                             merged_region_susie_alpha_res,
                                             merged_region_id_map){

  new_region_ids <- merged_region_id_map$region_id
  old_region_ids <- unlist(strsplit(merged_region_id_map$old_region_ids, ","))

  kept_finemap_res <- finemap_res[!finemap_res$region_id %in% old_region_ids, ]
  new_finemap_res <- merged_region_finemap_res[merged_region_finemap_res$region_id %in% new_region_ids, ]
  finemap_res <- rbind(kept_finemap_res, new_finemap_res)

  kept_susie_alpha_res <- susie_alpha_res[!susie_alpha_res$region_id %in% old_region_ids, ]
  new_susie_alpha_res <- merged_region_susie_alpha_res[merged_region_susie_alpha_res$region_id %in% new_region_ids, ]
  susie_alpha_res <- rbind(kept_susie_alpha_res, new_susie_alpha_res)

  finemap_res <- unique(finemap_res)
  rownames(finemap_res) <- NULL

  susie_alpha_res <- unique(susie_alpha_res)
  rownames(susie_alpha_res) <- NULL

  if (!setequal(finemap_res$id[finemap_res$group != "SNP"], susie_alpha_res$id))
    stop("finemap_res$id do not match with susie_alpha_res$id!")

  return(list("finemap_res" = finemap_res,
              "susie_alpha_res" = susie_alpha_res))
}


#' Updates cTWAS input data with merged region data
#'
#' @param region_data a list of original region data.
#' @param merged_region_data a list of merged region data.
#' @param region_info a data frame of original region definitions.
#' @param merged_region_info a data frame of original region definitions.
#' @param LD_map a data frame of original LD map.
#' @param merged_LD_map a data frame of merged LD map.
#' @param snp_map a list of original SNP info.
#' @param merged_snp_map a list of merged SNP info.
#' @param screened_region_L a vector of L for original screened regions.
#' @param merged_region_L a vector of L for merged regions.
#' @param merged_region_id_map a data frame of new region IDs and original regions IDs.
#'
#' @return a list with updated region_data, region_info, LD_map, snp_map, and L.
#'
#' @export
#'
update_merged_region_data <- function(region_data, merged_region_data,
                                      region_info, merged_region_info,
                                      LD_map, merged_LD_map,
                                      snp_map, merged_snp_map,
                                      screened_region_L, merged_region_L,
                                      merged_region_id_map){

  new_region_ids <- merged_region_id_map$region_id
  old_region_ids <- unlist(strsplit(merged_region_id_map$old_region_ids, ","))

  kept_region_info <- region_info[!region_info$region_id %in% old_region_ids, ]
  new_region_info <- merged_region_info[merged_region_info$region_id %in% new_region_ids, ]
  region_info <- rbind(kept_region_info, new_region_info)

  LD_map <- rbind(LD_map, merged_LD_map)
  idx <- match(region_info$region_id, LD_map$region_id)
  LD_map <- LD_map[idx, ]

  snp_map <- c(snp_map, merged_snp_map)
  snp_map <- snp_map[region_info$region_id]

  region_data <- c(region_data, merged_region_data)
  region_data <- region_data[region_info$region_id]

  kept_screened_region_L <- screened_region_L[!names(screened_region_L) %in% old_region_ids]
  new_region_L <- merged_region_L[new_region_ids]
  screened_region_L <- c(kept_screened_region_L, new_region_L)

  return(list("updated_region_data" = region_data,
              "updated_region_info" = region_info,
              "updated_LD_map" = LD_map,
              "updated_snp_map" = snp_map,
              "updated_region_L" = screened_region_L))
}

#' Updates cTWAS input data without LD with merged region data
#'
#' @param region_data a list of original region data.
#' @param merged_region_data a list of merged region data.
#' @param region_info a data frame of original region definitions.
#' @param merged_region_info a data frame of original region definitions.
#' @param snp_map a list of original SNP info.
#' @param merged_snp_map a list of merged SNP info.
#' @param merged_region_id_map a data frame of new region IDs and original regions IDs.
#'
#' @return a list with updated region_data, region_info, snp_map.
#'
#' @export
#'
update_merged_region_data_noLD <- function(region_data, merged_region_data,
                                           region_info, merged_region_info,
                                           snp_map, merged_snp_map,
                                           merged_region_id_map){

  new_region_ids <- merged_region_id_map$region_id
  old_region_ids <- unlist(strsplit(merged_region_id_map$old_region_ids, ","))

  kept_region_info <- region_info[!region_info$region_id %in% old_region_ids, ]
  new_region_info <- merged_region_info[merged_region_info$region_id %in% new_region_ids, ]
  region_info <- rbind(kept_region_info, new_region_info)

  snp_map <- c(snp_map, merged_snp_map)
  snp_map <- snp_map[region_info$region_id]

  region_data <- c(region_data, merged_region_data)
  region_data <- region_data[region_info$region_id]

  return(list("updated_region_data" = region_data,
              "updated_region_info" = region_info,
              "updated_snp_map" = snp_map))
}
