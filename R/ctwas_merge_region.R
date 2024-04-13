#' Identify overlapping regions
label_overlapping_regions <- function(boundary_genes) {
  # Sort the dataframe by 'start' and then by 'stop'
  boundary_genes <- boundary_genes[with(boundary_genes, order(chrom, region_start, region_stop)),]

  # Initialize the label counter and assign the first label
  current_label <- 1
  boundary_genes$merge_label <- current_label

  # Loop through the rows of the dataframe starting from the second row
  for (i in 2:nrow(boundary_genes)) {
    # Check if the current row overlaps with the previous row
    if (boundary_genes$region_start[i] <= boundary_genes$region_stop[i-1]) {
      boundary_genes$merge_label[i] <- current_label
    } else {
      current_label <- current_label + 1
      boundary_genes$merge_label[i] <- current_label
    }
  }

  return(boundary_genes)
}

#' #' Function to identify consecutive regions and label them
#' label_merged_regions <- function(df) {
#'   # Identify consecutive regions
#'   df$merged <- c(FALSE, df$start[-1] == df$stop[-length(df$stop)])
#'   df$label <- cumsum(!df$merged)
#'
#'   # Remove the consecutive column as it's no longer needed
#'   df$merged <- NULL
#'
#'   return(df)
#' }
#'
#'
#' merge_region <- function(regionlist, boundary_genes){
#'
#'   unique_regions <- unique(c(boundary_genes$region1,boundary_genes$region2))
#'   region_info <- data.frame(do.call(rbind,lapply(unique_regions, function(x){as.numeric(unlist(strsplit(x,"_")))})))
#'   colnames(region_info) <- c("chrom","start","stop")
#'   region_info <- region_info[order(region_info$chrom, region_info$start),]
#'
#'   regionlist_merged <- list()
#'   for(b in unique(region_info$chrom)){
#'     regions <- region_info[region_info$chrom==b,]
#'     regions <- label_merged_regions(regions)
#'
#'     for (i in unique(regions$label)){
#'       merged_region <- regions[regions$label==i,]
#'       region_id_new <- paste0(min(merged_region$start),"_",max(merged_region$stop))
#'       gid <- c()
#'       sid <- c()
#'       minpos <- c()
#'       maxpos <- c()
#'       LD_matrix <- c()
#'       SNP_info <- c()
#'       for(rn in 1:nrow(merged_region)){
#'         region_id <- paste0(merged_region$start[rn],"_",merged_region$stop[rn])
#'         region <- regionlist[[as.character(b)]][[region_id]]
#'         gid <- unique(c(gid, region[["gid"]]))
#'         sid <-  unique(c(sid, region[["sid"]]))
#'         minpos <- unique(c(minpos, region[["minpos"]]))
#'         maxpos <-  unique(c(maxpos, region[["maxpos"]]))
#'         LD_matrix <- unique(c(LD_matrix, region[["LD_matrix"]]))
#'         SNP_info <-  unique(c(SNP_info, region[["SNP_info"]]))
#'       }
#'       regionlist_merged[[as.character(b)]][[region_id_new]] <- list("gid"  = gid,
#'                                                 "sid"  = sid,
#'                                                 "start" = min(merged_region$start),
#'                                                 "stop" = max(merged_region$stop),
#'                                                 "minpos" = min(minpos),
#'                                                 "maxpos" = max(maxpos),
#'                                                 "LD_matrix" = LD_matrix,
#'                                                 "SNP_info" = SNP_info
#'                                                 )
#'     }
#'     loginfo("No. regions on chr%s after merging: %s", b, length(regionlist_merged[[as.character(b)]]))
#'   }
#'   return(regionlist_merged)
#' }
