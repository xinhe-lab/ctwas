
# converts old version regionlist to new version region_data
#' @importFrom logging loginfo
convert_regionlist_to_region_data <- function(regionlist,
                                              region_info,
                                              z_snp,
                                              z_gene,
                                              thin,
                                              ncore = 1){

  # convert the data structure from old version regionlist to new version region_data
  logging::loginfo("Convert regionlist to region_data")

  region_info <- region_info[order(region_info$chrom, region_info$start),]

  region_data <- list()
  for(b in 1:22){
    regionlist_chr <- regionlist[[b]]
    region_info_chr <- region_info[region_info$chrom == b, ]
    for(rn in names(regionlist_chr)){
      regiondata <- regionlist_chr[[rn]]
      regioninfo <- region_info_chr[as.numeric(rn),]
      if(regiondata$start != regioninfo$start || regiondata$stop != regioninfo$stop){
        stop(paste("Region does not match in chr:", b, "rn:", rn))
      }
      regiondata$chrom <- b
      region_id <- paste0(regiondata$chrom, "_", regiondata$start, "_", regiondata$stop)
      regiondata$region_id <- region_id
      regiondata$thin <- thin
      region_data[[region_id]] <- regiondata[c("region_id", "chrom", "start", "stop", "minpos", "maxpos", "thin", "gid", "sid")]
    }
  }
  logging::loginfo("%d regions in region_data", length(region_data))

  # add z-scores to region_data
  loginfo("Add region z-scores")
  region_data <- update_region_z(region_data, z_snp, z_gene, ncore = ncore)

  return(region_data)
}

# convert old version region tags to new version region ids
convert_region_tags_to_region_id <- function(region_info, region_tags1, region_tags2){

  stopifnot(length(region_tags1) == length(region_tags2))
  region_tags <- paste0(region_tags1, "_", region_tags2)
  region_tags_df <- unique(data.frame(region_tag1 = region_tags1,
                                      region_tag2 = region_tags2,
                                      region_tag = region_tags))

  region_info <- region_info[order(region_info$chrom, region_info$start),]
  if (is.null(region_info$region_id)) {
    region_info$region_id <- paste0(region_info$chrom, "_", region_info$start, "_", region_info$stop)
  }

  region_ids <- sapply(1:nrow(region_tags_df), function(i){
    b <- as.integer(region_tags_df$region_tag1[i])
    rn <- as.integer(region_tags_df$region_tag2[i])
    region_info_chr <- region_info[region_info$chrom == b, ]
    region_id <- region_info_chr[rn, "region_id"]
    region_id
  })
  region_tags_df$region_id <- region_ids
  idx <- match(region_tags, region_tags_df$region_tag)
  region_ids <- region_tags_df$region_id[idx]
  return(region_ids)
}
