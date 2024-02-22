
#' filter regions based on probability of at most 1 causal effect
filter_regions <- function(regionlist, group_prior, prob_single = 0.8, zdf){
  regionlist2 <- regionlist
  for (b in 1: length(regionlist)){
    for (rn in names(regionlist[[b]])) {
      gid <- regionlist[[b]][[rn]][["gid"]]
      sid <- regionlist[[b]][[rn]][["sid"]]
      gs_type <- zdf$type[match(c(gid,sid), zdf$id)]

      group_size <- table(gs_type)[names(group_prior)]
      group_size[is.na(group_size)] <- 0

      P1 <- prod((1-group_prior)^group_size) * (1 + sum(group_size*(group_prior/(1-group_prior))))

      if (P1 < prob_single){
        regionlist2[[b]][[rn]] <- NULL
      }
    }
  }
  regionlist2
}

#' assign regions to cores
#'
#' @param ncore integer, numeber of cores, at least 1
#' regions allocated to given number of cores
#' regionlist need to contain at least 1 non-empty
#'
region2core <- function(regionlist, ncore = 1){
  dflist <- list()
  for (b in 1:length(regionlist)){
    if (length(regionlist[[b]]) > 0){
      dflist[[b]] <- data.frame("b" = b, "rn"= names(regionlist[[b]]), stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, dflist)
  if (ncore > 1) {
    d <- cut(1:nrow(df), ncore, labels = FALSE)
    corelist <- split(df,d)
  } else {
    corelist <- list("1" = df)
  }
  return(corelist)
}

# select and assemble a subset of regionlist by region_tags
subset_regionlist <- function(regionlist, region_tags = NULL){

  if (length(region_tags) == 0){
    stop("No regions selected!")
  }

  # subset regionlist by region_tags
  region_tags_df <- data.frame(
    b = sapply(strsplit(region_tags, "_"), "[[", 1),
    rn = sapply(strsplit(region_tags, "_"), "[[", 2))

  regionlist_subset <- vector("list", length = 22)
  for (b in region_tags_df$b) {
    rn <- region_tags_df[region_tags_df$b == b, "rn"]
    regionlist_subset[[b]] <- regionlist[[b]][as.character(rn)]
  }

  # coordinates of the subset regions
  temp_regs <- lapply(1:22, function(x) cbind(x,
                                              unlist(lapply(regionlist_subset[[x]], "[[", "start")),
                                              unlist(lapply(regionlist_subset[[x]], "[[", "stop"))))

  regs_subset <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
  loginfo("Subset %d regions from the regionlist", nrow(regs_subset))

  return(list("regionlist" = regionlist_subset,
              "regs" = regs))
}


#' Merge regions.
#'
#' @param region_info a data frame of region definition and associated file names
#'
#' @param regionlist the list of regions to be merged.
#'
#' @param chrom the list of chromosomes to be merged.
#'
#' @importFrom logging loginfo
#'
#' @return a list of merged regions.
#'
merge_regions <- function(region_info, regionlist, chrom = 1:22) {

  for (b in chrom){

    for (rn in 2:nrow(region_info)){
      current <- regionlist[[b]][[as.character(rn)]]
      previous <- regionlist[[b]][[as.character(rn - 1)]]

      if (length(intersect(current[["gid"]], previous[["gid"]]))> 0){

        gidx <-  unique(c(previous[["gidx"]], current[["gidx"]]))
        sidx <-  unique(c(previous[["sidx"]], current[["sidx"]]))

        gid <- geneinfo$id[gidx]
        sid <- snpinfo$id[sidx]

        regionlist[[b]][[as.character(rn)]] <- list("gidx" = gidx,
                                                    "gid"  = gid,
                                                    "sidx" = sidx,
                                                    "sid"  = sid,
                                                    "start" = previous$start,
                                                    "stop" = current$stop,
                                                    "minpos" = previous$minpos,
                                                    "maxpos" = current$maxpos)

        regionlist[[b]][[as.character(rn -1)]] <- NULL
      }
    }
  }

  loginfo("No. regions with at least one SNP/gene for chr%s after merging: %s",
          b, length(regionlist[[b]]))

  ld_Rf <- ld_Rfs[b]
  ld_Rinfo <- as.data.frame(data.table::fread(ld_Rf, header = T))

  if (!isTRUE(merge) & nrow(regions) >=2){
    for (rn in 1:(nrow(regions)-1)){
      current <- regionlist[[b]][[as.character(rn)]]
      nextone <- regionlist[[b]][[as.character(rn+1)]]
      gnames <- regionlist[[b]][[as.character(rn)]][["gid"]]
      tmp_region <- regionlist[[b]][[as.character(rn)]]
      if(length(gnames>0)){
        ifreg <- ifelse(regionlist[[b]][[as.character(rn)]][["start"]] < ld_Rinfo[, "stop"] & regionlist[[b]][[as.character(rn)]][["stop"]] >= ld_Rinfo[, "start"], T, F)
        regRDS <- ld_Rinfo[ifreg, "RDS_file"]
        R_snp_anno <- as.data.frame(do.call(rbind, lapply(regRDS, ctwas:::read_ld_Rvar_RDS)))
        for (i in 1:length(gnames)){
          gname <- gnames[i]
          wgt <- wgtlistall[[gname]]
          snpnames <- rownames(wgt)
          ld.idx <- match(snpnames, R_snp_anno$id)
          if(anyNA(ld.idx)){
            thisindex <- !is.na(ld.idx)
            nextindex <- is.na(ld.idx)
            thisr2 <- sum(wgt[thisindex]^2)
            nextr2 <- sum(wgt[nextindex]^2)
            if(thisr2<nextr2){
              #modify weights file - drop weights in other regions
              tmp_wgt <- wgtlistall[[gname]][nextindex]
              if(length(tmp_wgt)==1){
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[nextindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
              else{
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[nextindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
              #add gene to next region
              regionlist[[b]][[as.character(rn+1)]][["gidx"]] <- c(regionlist[[b]][[as.character(rn+1)]][["gidx"]],tmp_region[["gidx"]][which(gnames==gname)])
              regionlist[[b]][[as.character(rn+1)]][["gid"]] <- c(regionlist[[b]][[as.character(rn+1)]][["gid"]],gname)
              #remove gene from this region
              regionlist[[b]][[as.character(rn)]][["gidx"]] <- regionlist[[b]][[as.character(rn)]][["gidx"]][which(gnames!=gname)]
              regionlist[[b]][[as.character(rn)]][["gid"]] <- regionlist[[b]][[as.character(rn)]][["gid"]][!regionlist[[b]][[as.character(rn)]][["gid"]]==gname]
            }
            else{
              #modify weights file - drop weights in other regions
              tmp_wgt <- wgtlistall[[gname]][thisindex]
              if(length(tmp_wgt)==1){
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = 1,ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[thisindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
              else{
                wgtlistall[[gname]] <- matrix(tmp_wgt,nrow = length(tmp_wgt),ncol = 1)
                rownames(wgtlistall[[gname]]) <- snpnames[thisindex]
                colnames(wgtlistall[[gname]]) <- "weight"
              }
            }
          }
        }
      }
    }
  }

  return(regionlist)
}
