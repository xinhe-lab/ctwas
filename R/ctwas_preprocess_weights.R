
#' Preprocess PredictDB/FUSION weights and harmonize with LD reference
#'
#' @param weight_file a string or a vector, pointing path to one or multiple sets of weights in PredictDB or FUSION format.
#'
#' @param region_info a data frame of region definition and associated file names.
#'
#' @param gwas_snp_ids a vector of SNP IDs in GWAS summary statistics (z_snp$id).
#'
#' @param snp_info a list of SNP info data frames for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref", and "region_id".
#'
#' @param type a string, specifying QTL type of each weight file, e.g. expression, splicing, protein.
#'
#' @param context a string, specifying tissue/cell type/condition of each weight file, e.g. Liver, Lung, Brain.
#'
#' @param weight_format a string, specifying format of each weight file, e.g. PredictDB, FUSION.
#'
#' @param filter_protein_coding_genes TRUE/FALSE. If TRUE, keep protein coding genes only. This option is only for PredictDB weights
#'
#' @param load_predictdb_LD TRUE/FALSE. If TRUE, load pre-computed LD among weight SNPs. This option is only for PredictDB weights
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @param method_FUSION a string, specifying the method to choose in FUSION models
#'
#' @param genome_version a string, specifying the genome version of FUSION models
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom foreach %dopar% foreach
#'
#' @return a list of processed weights
#'
#' @export
#'
preprocess_weights <- function(weight_file,
                               region_info,
                               gwas_snp_ids,
                               snp_info,
                               type,
                               context,
                               weight_format = c("PredictDB", "FUSION"),
                               ncore = 1,
                               drop_strand_ambig = TRUE,
                               scale_by_ld_variance = FALSE,
                               filter_protein_coding_genes = FALSE,
                               load_predictdb_LD = FALSE,
                               method_FUSION = c("lasso","enet","top1","blup"),
                               genome_version = c("b38","b37"),
                               logfile = NULL){

if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }

  # check input arguments
  weight_format <- match.arg(weight_format)
  method_FUSION <- match.arg(method_FUSION)
  genome_version <- match.arg(genome_version)

  loginfo("Load weight: %s", weight_file)

  if (length(weight_file) > 1) {
    stop("Please provide only one weight file in weight_file.")
  }
  stopifnot(file.exists(weight_file))

  # Check LD reference SNP info
  if (class(snp_info) == "list") {
    snp_info <- as.data.frame(data.table::rbindlist(snp_info, idcol = "region_id"))
  }
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info))){
    stop("SNP info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  # set default type and context
  if (missing(type)) {
    type <- "gene"
  }

  if (missing(context)) {
    context <- tools::file_path_sans_ext(basename(weight_file))
  }

  loginfo("type: %s", type)
  loginfo("context: %s", context)

  weights <- list()
  loaded_weight <- load_weights(weight_file,
                                weight_format,
                                filter_protein_coding_genes = filter_protein_coding_genes,
                                load_predictdb_LD = load_predictdb_LD,
                                method_FUSION = method_FUSION,
                                genome_version = genome_version,
                                ncore=ncore)
  weight_table <- loaded_weight$weight_table
  weight_name <- loaded_weight$weight_name
  R_wgt_all <- loaded_weight$R_wgt
  if(!is.null(R_wgt_all)){
    weight_table <- weight_table[weight_table$gene %in% unique(R_wgt_all$GENE), ] #remove genes without predictdb LD
  }
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %d in %s", length(gnames), weight_name)
  # remove variants in weight table, but not in LD reference and GWAS
  loginfo("Number of variants in weights: %d", length(unique(weight_table$rsid)))
  # take the intersect of SNPs in weights, LD reference and SNPs in z_snp
  snpnames <- Reduce(intersect, list(weight_table$rsid, snp_info$id, gwas_snp_ids))
  # loginfo("Remove %d variants after intersecting with LD reference and GWAS", length(setdiff(weight_table$rsid, snpnames)))
  weight_table <- weight_table[weight_table$rsid %in% snpnames, ]
  # loginfo("Remove %s genes after intersecting with LD reference and GWAS", length(setdiff(gnames, weight_table$gene)))
  gnames <- unique(weight_table$gene)
  loginfo("%d variants and %d genes left after intersecting with LD reference and z_snp", length(snpnames), length(gnames))
  # subset to variants in weight table
  snp_info_wgt <- snp_info[snp_info$id %in% weight_table$rsid,]
  loginfo("Harmonizing weights with LD reference ...")
  rsid_varID <- weight_table[,c("rsid", "varID")]

  for (i in 1:length(gnames)){
    gname <- gnames[i]
    wgt <- weight_table[weight_table$gene==gname,]
    wgt.matrix <- as.matrix(wgt[, "weight", drop = F])
    rownames(wgt.matrix) <- wgt$rsid
    chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))
    chrom <- unique(as.integer(gsub("chr", "", chrpos[, 1])))
    if (length(chrom) > 1) {
      stop("More than one chrom in weight for %s!", gname)
    }
    snp_pos <- as.integer(chrpos[, 2])
    wgt_ld_idx <- match(wgt$rsid, snp_info_wgt$id)
    snp_info_pos <- as.integer(snp_info_wgt$pos[wgt_ld_idx])
    if (any(snp_pos != snp_info_pos)){
      warning(sprintf("Variant positions in %s weights are different from positions in snp_info. Use the positions in snp_info instead.", gname))
      snp_pos <- snp_info_pos
    }
    snps <- data.frame(chrom = chrom,
                       id = wgt$rsid,
                       cm = "0",
                       pos = snp_pos,
                       alt = wgt$eff_allele,
                       ref = wgt$ref_allele,
                       stringsAsFactors = F)

    w <- harmonize_weights(wgt.matrix,
                          snps,
                          snp_info_wgt,
                          drop_strand_ambig = drop_strand_ambig)
    wgt.matrix <- w[["wgt"]]
    snps <- w[["snps"]]
    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, "weight"]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]

    snpnames <- intersect(rownames(wgt.matrix), snp_info_wgt$id)
    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    wgt <- wgt.matrix[wgt.idx, "weight", drop = F]

    snps.idx <- match(snpnames, snps$id)
    snps <- snps[snps.idx,]

    if (scale_by_ld_variance){
      snp_info_wgt.idx <- match(snpnames, snp_info_wgt$id)
      wgt <- wgt*sqrt(snp_info_wgt$variance[snp_info_wgt.idx])
    }

    n_wgt <- nrow(wgt.matrix)

    if(n_wgt>0){
      p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
      p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
      # weight_id <- paste0(gname, "|", type, "|", context)
      weight_id <- paste0(gname, "|", weight_name)

      #Add LD matrix of weights
      if(!is.null(R_wgt_all)){
        R_wgt <- get_weight_LD(R_wgt_all,gname,rsid_varID)
        R_wgt <- R_wgt[snps$id, snps$id, drop=F]
      }
      else{
        R_wgt <- NULL
      }
      weights[[weight_id]] <- list(chrom=chrom, p0=p0, p1=p1, wgt=wgt, R_wgt=R_wgt,
                                   gene_name=gname, weight_name=weight_name,
                                   type = type, context = context, n_wgt=n_wgt)
    }
  }

  if(!load_predictdb_LD){
    loginfo("Computing LD between variants in weights ...")
    weight_info <- as.data.frame(do.call(rbind, weights)[,c("chrom","p0","p1","gene_name","weight_name","type","context")])
    weight_info$weight_id <- paste0(weight_info$gene_name, "|", weight_name)
    for (k in 1:nrow(weight_info)) {
      chrom <- weight_info[k, "chrom"]
      p0 <- weight_info[k, "p0"]
      p1 <- weight_info[k, "p1"]
      idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
      weight_info[k, "region_id"] <- paste(sort(region_info[idx, "region_id"]), collapse = ";")
    }
    # impute LD for weights for each chromosome
    cl <- parallel::makeCluster(ncore, outfile = "")
    doParallel::registerDoParallel(cl)

    for (b in unique(weight_info$chrom)) {
      loginfo("Computing LD for weight variants on chr%s", b)
      weightinfo <- weight_info[weight_info$chrom == b, ]

      if (nrow(weightinfo) > 0) {
        batches <- names(sort(-table(weightinfo$region_id)))
        corelist <- lapply(1:ncore, function(core){
          batches_core <- batches[0:ceiling(length(batches)/ncore-1)*ncore+core];
          batches_core[!is.na(batches_core)]})
        names(corelist) <- 1:ncore

        outlist <- foreach(core = 1:ncore, .combine = "c", .packages = c("ctwas")) %dopar% {
          batches <- corelist[[core]]
          outlist_core <- list()

          for (batch in batches) {
            # load the R_snp and SNP info for the region
            region_ids <- strsplit(batch, ";")[[1]]
            reg_idx <- match(region_ids, region_info$region_id)
            if (length(reg_idx) > 1){
              R_snp <- lapply(region_info$LD_matrix[reg_idx], load_LD)
              R_snp <- suppressWarnings({as.matrix(Matrix::bdiag(R_snp))})
            }
            else{
              R_snp <- load_LD(region_info$LD_matrix[reg_idx])
            }
            snpinfo <- read_snp_info_files(region_info$SNP_info[reg_idx])
            weight_ids <- weightinfo[weightinfo$region_id == batch, "weight_id"]
            for(weight_id in weight_ids){
              snpnames <- rownames(weights[[weight_id]]$wgt)
              sidx <- match(snpnames, snpinfo$id)
              R_wgt <- R_snp[sidx, sidx, drop=F]
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
