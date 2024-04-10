
#' Preprocess PredictDB/Fusion weights and harmonize with LD reference
#'
#' @param weight_file a string or a vector, pointing path to one or multiple sets of weights in PredictDB or Fusion format.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param z_snp A data frame with columns: "id", "A1", "A2", "z". giving the z scores for
#' snps. "A1" is effect allele. "A2" is the other allele. For harmonized data, A1 and A2 are not required.
#'
#' @param type a string or a vector, specifying QTL type of each weight file, e.g. eQTL, sQTL, pQTL.
#'
#' @param context a string or a vector, specifying tissue/cell type/condition of each weight file, e.g. Liver, Lung, Brain.
#'
#' @param weight_format a string or a vector, specifying format of each weight file, e.g. PredictDB, Fusion.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD among weight SNPs. This option is only for PredictDB weights
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param method_Fusion a string, specifying the method to choose in Fusion models
#'
#' @return a list of processed weights
#'
#' @export
#'
preprocess_weights <- function(weight_file,
                               region_info,
                               z_snp,
                               type = NULL,
                               context = NULL,
                               weight_format = c("PredictDB", "Fusion"),
                               ncore = 1,
                               drop_strand_ambig = TRUE,
                               filter_protein_coding_genes = FALSE,
                               load_predictdb_LD = FALSE,
                               method_Fusion = "enet"){

  # check input arguments
  weight_format <- match.arg(weight_format)
  if (is.null(type)) {
    type <- "gene"
  }
  if (is.null(context)) {
    context <- tools::file_path_sans_ext(basename(weight_file))
  }

  loginfo("Load weight: %s", weight_file)
  loginfo("type: %s", type)
  loginfo("context: %s", context)

  # load LD SNPs information
  region_info <- region_info[order(region_info$chrom, region_info$start),]
  ld_snpinfo <- read_LD_SNP_files(region_info$SNP_info)

  weights <- list()
  loaded_weight <- load_weights(weight_file,
                                weight_format,
                                ld_snpinfo,
                                filter_protein_coding_genes = filter_protein_coding_genes,
                                load_predictdb_LD = load_predictdb_LD,
                                method_Fusion = method_Fusion,
                                ncore=ncore)
  weight_table <- loaded_weight$weight_table
  weight_name <- loaded_weight$weight_name
  R_wgt_all <- loaded_weight$R_wgt
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %d in %s", length(gnames), weight_name)
  # remove variants in weight table, but not in LD reference and GWAS
  loginfo("Number of variants in weights: %d", length(unique(weight_table$rsid)))
  # take the intersect of SNPs in weights, LD reference and z_snp
  snpnames <- Reduce(intersect, list(weight_table$rsid, ld_snpinfo$id, z_snp$id))
  # loginfo("Remove %d variants after intersecting with LD reference and GWAS", length(setdiff(weight_table$rsid, snpnames)))
  weight_table <- weight_table[weight_table$rsid %in% snpnames, ]
  # loginfo("Remove %s genes after intersecting with LD reference and GWAS", length(setdiff(gnames, weight_table$gene)))
  gnames <- unique(weight_table$gene)
  loginfo("%d variants and %d genes left after intersecting with LD reference and GWAS z_snp", length(snpnames), length(gnames))
  # subset to variants in weight table
  ld_snpinfo_wgt <- ld_snpinfo[ld_snpinfo$id %in% weight_table$rsid,]
  loginfo("Harmonize weights with LD reference")
  rsid_varID <- weight_table[,c("rsid", "varID")]

  pb <- txtProgressBar(min = 0, max = length(gnames), initial = 0, style = 3)
  for (j in 1:length(gnames)){
    gname <- gnames[j]
    wgt <- weight_table[weight_table$gene==gname,]
    wgt.matrix <- as.matrix(wgt[, "weight", drop = F])
    rownames(wgt.matrix) <- wgt$rsid
    chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))
    chrom <- unique(as.integer(gsub("chr", "", chrpos[, 1])))
    if (length(chrom) > 1) {
      stop("More than one chrom in weight for %s!", gname)
    }
    wgt_ld_idx <- match(wgt$rsid, ld_snpinfo_wgt$id)
    snp_pos <- as.integer(chrpos[, 2])
    ld_snp_pos <- as.integer(ld_snpinfo_wgt[wgt_ld_idx, "pos"])
    if (!all.equal(snp_pos, ld_snp_pos)){
      loginfo("Variant positions in weights are not the same as positions in LD reference. Use positions in LD reference.")
      snp_pos <- ld_snp_pos
    }
    snps <- data.frame(chrom = chrom,
                       id = wgt$rsid,
                       cm = "0",
                       pos = snp_pos,
                       alt = wgt$eff_allele,
                       ref = wgt$ref_allele,
                       stringsAsFactors = F)

    w <- harmonize_wgt_ld(wgt.matrix,
                          snps,
                          ld_snpinfo_wgt,
                          drop_strand_ambig = drop_strand_ambig)
    wgt.matrix <- w[["wgt"]]
    snps <- w[["snps"]]
    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, "weight"]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]

    snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo_wgt$id)
    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    wgt <- wgt.matrix[wgt.idx, "weight", drop = F]

    snps.idx <- match(snpnames, snps$id)
    snps <- snps[snps.idx,]

    #if (isTRUE(scale_by_ld_variance)){
    #  ld_snpinfo_wgt.idx <- match(snpnames, ld_snpinfo_wgt$id)
    #  wgt <- wgt*sqrt(ld_snpinfo_wgt$variance[ld_snpinfo_wgt.idx])
    #}

    n_wgt <- nrow(wgt.matrix)

    if(n_wgt>0){
      p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
      p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
      # weight_id <- paste0(gname, "|", type, "|", context)
      weight_id <- paste0(gname, "|", weight_name)

      #Add LD matrix of weights
      if(!is.null(R_wgt_all)){
        R_wgt <- get_weight_LD(R_wgt_all,gname,rsid_varID)
        R_wgt <- R_wgt[snps$id,snps$id,drop=F]
      }
      else{
        R_wgt <- NULL
      }
      weights[[weight_id]] <- list(chrom=chrom, p0=p0, p1=p1, wgt=wgt, R_wgt=R_wgt,
                                   gene_name=gname, weight_name=weight_name,
                                   type = type, context = context, n_wgt=n_wgt)
    }
    setTxtProgressBar(pb, j)
  }
  close(pb)

  if(!load_predictdb_LD){
    loginfo("Computing LD between variants in weights.")
    weight_info <- as.data.frame(do.call(rbind, weights)[,c("chrom","p0","p1","gene_name","weight_name","type","context")])
    # weight_info$weight_id <- paste0(weight_info$gene_name, "|", weight_info$type, "|", weight_info$context)
    weight_info$weight_id <- paste0(weight_info$gene_name, "|", weight_name)
    for (k in 1:nrow(weight_info)) {
      chrom <- weight_info[k, "chrom"]
      p0 <- weight_info[k, "p0"]
      p1 <- weight_info[k, "p1"]
      idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
      weight_info[k, "region_tag"] <- paste(sort(region_info[idx, "region_tag"]), collapse = ";")
    }
    # impute LD for weights for each chromosome
    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)

    for (b in unique(weight_info$chrom)) {
      loginfo("Computing LD for weight variants on chr%s", b)
      weightinfo <- weight_info[weight_info$chrom == b, ]

      if (nrow(weightinfo) > 0) {
        batches <- names(sort(-table(weightinfo$region_tag)))
        corelist <- lapply(1:ncore, function(core){
          batches_core <- batches[0:ceiling(length(batches)/ncore-1)*ncore+core];
          batches_core[!is.na(batches_core)]})
        names(corelist) <- 1:ncore

        outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas")) %dopar% {
          batches <- corelist[[core]]
          outlist_core <- list()

          for (batch in batches) {
          # load the R_snp and SNP info for the region
            region_tags <- strsplit(batch, ";")[[1]]
            reg_idx <- match(region_tags, region_info$region_tag)
            if (length(reg_idx) > 1){
              R_snp <- lapply(region_info$LD_matrix[reg_idx], load_LD)
              R_snp <- suppressWarnings({as.matrix(Matrix::bdiag(R_snp))})
            }
            else{
              R_snp <- load_LD(region_info$LD_matrix[reg_idx])
            }
            R_snpinfo <- read_LD_SNP_files(region_info$SNP_info[reg_idx])
            weight_ids <- weightinfo[weightinfo$region_tag == batch, "weight_id"]
            for(weight_id in weight_ids){
              snpnames <- rownames(weights[[weight_id]]$wgt)
              ld.idx <- match(snpnames, R_snpinfo$id)
              R_wgt <- R_snp[ld.idx, ld.idx, drop=F]
              rownames(R_wgt) <- snpnames
              colnames(R_wgt) <- snpnames
              outlist_core[[weight_id]] <- R_wgt
            }
          }
          outlist_core
        }
        for(l in names(outlist)){
          weights[[l]][["R_wgt"]] <- outlist[[l]]
        }
      }
    }
  }
  return(weights)
}
