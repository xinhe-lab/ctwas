merge_region <- function(regionlist, boundary_genes){

  unique_regions <- unique(c(boundary_genes$region1,boundary_genes$region2))
  region_info <- data.frame(do.call(rbind,lapply(unique_regions, function(x){as.numeric(unlist(strsplit(x,"_")))})))
  colnames(region_info) <- c("chrom","start","stop")
  region_info <- region_info[order(region_info$chrom, region_info$start),]

  regionlist_merged <- list()
  for(b in unique(region_info$chrom)){
    regions <- region_info[region_info$chrom==b,]
    regions <- label_merged_regions(regions)

    for (i in unique(regions$label)){
      merged_region <- regions[regions$label==i,]
      region_id_new <- paste0(min(merged_region$start),"_",max(merged_region$stop))
      gid <- c()
      sid <- c()
      minpos <- c()
      maxpos <- c()
      LD_matrix <- c()
      SNP_info <- c()
      for(rn in 1:nrow(merged_region)){
        region_id <- paste0(merged_region$start[rn],"_",merged_region$stop[rn])
        region <- regionlist[[as.character(b)]][[region_id]]
        gid <- unique(c(gid, region[["gid"]]))
        sid <-  unique(c(sid, region[["sid"]]))
        minpos <- unique(c(minpos, region[["minpos"]]))
        maxpos <-  unique(c(maxpos, region[["maxpos"]]))
        LD_matrix <- unique(c(LD_matrix, region[["LD_matrix"]]))
        SNP_info <-  unique(c(SNP_info, region[["SNP_info"]]))
      }
      regionlist_merged[[as.character(b)]][[region_id_new]] <- list("gid"  = gid,
                                                "sid"  = sid,
                                                "start" = min(merged_region$start),
                                                "stop" = max(merged_region$stop),
                                                "minpos" = min(minpos),
                                                "maxpos" = max(maxpos),
                                                "LD_matrix" = LD_matrix,
                                                "SNP_info" = SNP_info
                                                )
    }
    loginfo("No. regions on chr%s after merging: %s", b, length(regionlist_merged[[as.character(b)]]))
  }
  return(regionlist_merged)
}
