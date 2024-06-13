
#' @title Preprocess PredictDB/FUSION weights and harmonize with LD reference
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
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{load_predictdb_LD = FALSE}.
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
#' @param fusion_genome_version a string, specifying the genome version of FUSION models
#'
#' @param fusion_top_n_snps a number, specifying the top n weight SNPs included in FUSION models. If NULL, using all weight SNPs
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader()} function to load LD matrix.
#'
#' @param LD_loader a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param logfile the log file, if NULL will print log info on screen.
#'
#' @return a list of processed weights
#'
#' @importFrom utils head
#' @importFrom logging addHandler loginfo writeToFile
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel mclapply makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table rbindlist
#' @importFrom tools file_path_sans_ext
#' @importFrom stats complete.cases
#' @importFrom Matrix bdiag
#'
#' @export
#'
preprocess_weights <- function(weight_file,
                               region_info,
                               gwas_snp_ids,
                               type,
                               context,
                               snp_info,
                               LD_info = NULL,
                               weight_format = c("PredictDB", "FUSION"),
                               ncore = 1,
                               drop_strand_ambig = TRUE,
                               scale_by_ld_variance = FALSE,
                               filter_protein_coding_genes = FALSE,
                               load_predictdb_LD = FALSE,
                               method_FUSION = c("lasso","enet","top1","blup"),
                               fusion_genome_version = c("b38","b37"),
                               fusion_top_n_snps = NULL,
                               LD_format = c("rds", "rdata", "csv", "txt", "custom"),
                               LD_loader = NULL,
                               logfile = NULL){
  if (!is.null(logfile)) {
    addHandler(writeToFile, file = logfile, level = "DEBUG")
  }
  # check input arguments
  weight_format <- match.arg(weight_format)
  method_FUSION <- match.arg(method_FUSION)
  fusion_genome_version <- match.arg(fusion_genome_version)
  LD_format <- match.arg(LD_format)

  if (length(weight_file) > 1) {
    stop("Please provide only one weight file in weight_file.")
  }
  stopifnot(file.exists(weight_file))

  # Check LD reference SNP info
  if (inherits(snp_info,"list")) {
    snp_info_df <- as.data.frame(rbindlist(snp_info, idcol = "region_id"))
  }
  target_header <- c("chrom", "id", "pos", "alt", "ref")
  if (!all(target_header %in% colnames(snp_info_df))){
    stop("SNP info needs to contain the following columns: ",
         paste(target_header, collapse = " "))
  }

  # set default type and context
  if (missing(type)) {
    type <- "gene"
  }

  if (missing(context)) {
    context <- file_path_sans_ext(basename(weight_file))
  }

  loginfo("Load weight: %s", weight_file)
  loginfo("type: %s", type)
  loginfo("context: %s", context)

  weights <- list()
  loaded_weight <- load_weights(weight_file,
                                weight_format,
                                filter_protein_coding_genes = filter_protein_coding_genes,
                                load_predictdb_LD = load_predictdb_LD,
                                method_FUSION = method_FUSION,
                                fusion_genome_version = fusion_genome_version,
                                ncore=ncore)

  weight_table <- loaded_weight$weight_table
  weight_name <- loaded_weight$weight_name
  R_wgt_all <- loaded_weight$R_wgt
  if (!is.null(R_wgt_all)) {
    weight_table <- weight_table[weight_table$gene %in% unique(R_wgt_all$GENE), ] #remove genes without predictdb LD
  }
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %d in %s", length(gnames), weight_name)
  # remove variants in weight table, but not in LD reference and GWAS
  loginfo("Number of variants in weights: %d", length(unique(weight_table$rsid)))
  # take the intersect of SNPs in weights, LD reference and SNPs in z_snp
  snpnames <- Reduce(intersect, list(weight_table$rsid, snp_info_df$id, gwas_snp_ids))
  # loginfo("Remove %d variants after intersecting with LD reference and GWAS", length(setdiff(weight_table$rsid, snpnames)))
  weight_table <- weight_table[weight_table$rsid %in% snpnames, ]
  # loginfo("Remove %s genes after intersecting with LD reference and GWAS", length(setdiff(gnames, weight_table$gene)))
  gnames <- unique(weight_table$gene)
  loginfo("%d variants and %d genes left after filtering by GWAS and reference SNPs", length(snpnames), length(gnames))
  # subset to variants in weight table
  snp_info_wgt <- snp_info_df[snp_info_df$id %in% weight_table$rsid,]
  loginfo("Harmonizing weights with LD reference ...")
  rsid_varID <- weight_table[,c("rsid", "varID")]

  for (i in 1:length(gnames)) {
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
    if (any(snp_pos != snp_info_pos)) {
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

    if (weight_format == "FUSION") {
      wgt <- wgt[order(-abs(wgt[,"weight"])), , drop = F]
      if (!is.null(fusion_top_n_snps)) {
        wgt <- head(wgt,fusion_top_n_snps)
      }
      snpnames <- intersect(rownames(wgt), snp_info_wgt$id)
    }

    snps.idx <- match(snpnames, snps$id)
    snps <- snps[snps.idx,]

    if (scale_by_ld_variance){
      snp_info_wgt.idx <- match(snpnames, snp_info_wgt$id)
      wgt <- wgt*sqrt(snp_info_wgt$variance[snp_info_wgt.idx])
    }

    n_wgt <- nrow(wgt)

    if (n_wgt>0) {
      p0 <- min(snps[snps[, "id"] %in% snpnames, "pos"])
      p1 <- max(snps[snps[, "id"] %in% snpnames, "pos"])
      # weight_id <- paste0(gname, "|", type, "|", context)
      weight_id <- paste0(gname, "|", weight_name)

      # Add LD matrix of weights
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

  if (!load_predictdb_LD) {
    loginfo("Computing LD between variants in weights ...")
    weights <- compute_weight_LD_from_ref(weights,
                                          weight_name,
                                          region_info = region_info,
                                          LD_info = LD_info,
                                          snp_info = snp_info,
                                          LD_format = LD_format,
                                          LD_loader=LD_loader,
                                          ncore = ncore)
  }

  return(weights)
}


# gets pre-computed LD matrix from predictedDB weights
#' @importFrom stats setNames
get_weight_LD <- function (R_wgt_all, gname, rsid_varID){
  R_wgt <- R_wgt_all[R_wgt_all$GENE == gname,]
  #convert covariance to correlation
  R_wgt_stdev <- R_wgt[R_wgt$RSID1==R_wgt$RSID2,]
  R_wgt_stdev <- setNames(sqrt(R_wgt_stdev$VALUE), R_wgt_stdev$RSID1)
  R_wgt$VALUE <- R_wgt$VALUE/(R_wgt_stdev[R_wgt$RSID1]*R_wgt_stdev[R_wgt$RSID2])

  unique_id <- unique(c(R_wgt$RSID1, R_wgt$RSID2))

  # Create an empty correlation matrix
  n <- length(unique_id)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)

  # Fill in the correlation values
  for (i in 1:n) {
    for (j in i:n) {  # Only iterate over half of the matrix
      if (i == j) {
        cor_matrix[i, j] <- 1  # Diagonal elements are 1
      } else {
        # Check if there are any matches for the RSID combination
        matches <- R_wgt[R_wgt$RSID1 == unique_id[i] & R_wgt$RSID2 == unique_id[j], "VALUE"]
        if (length(matches) > 0) {
          cor_matrix[i, j] <- matches
          cor_matrix[j, i] <- matches  # Set symmetric value
        } else {
          cor_matrix[i, j] <- NA  # No correlation value found
          cor_matrix[j, i] <- NA  # No correlation value found
        }
      }
    }
  }

  rownames(cor_matrix) <- rsid_varID$rsid[match(unique_id, rsid_varID$varID)]
  colnames(cor_matrix) <- rsid_varID$rsid[match(unique_id, rsid_varID$varID)]

  return(cor_matrix)
}

# compute LD for weight variants using reference LD
#' @importFrom parallel mclapply
#' @importFrom Matrix bdiag
#' @importFrom logging loginfo
compute_weight_LD_from_ref <- function(weights,
                                       weight_name,
                                       region_info,
                                       LD_info,
                                       snp_info,
                                       LD_format = c("rds", "rdata", "csv", "txt", "custom"),
                                       LD_loader = NULL,
                                       ncore = 1) {

  if (is.null(LD_info) || is.null(snp_info)) {
    stop("LD_info and snp_info are required for computing LD")
  }

  LD_format <- match.arg(LD_format)

  weight_info <- lapply(names(weights), function(x){
    as.data.frame(weights[[x]][c("chrom", "p0","p1", "gene_name", "weight_name", "type","context")])})
  weight_info <- do.call(rbind, weight_info)
  weight_info$weight_id <- paste0(weight_info$gene_name, "|", weight_name)
  # get the regions overlapping with each gene
  for (k in 1:nrow(weight_info)) {
    chrom <- weight_info[k, "chrom"]
    p0 <- weight_info[k, "p0"]
    p1 <- weight_info[k, "p1"]
    idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
    weight_info[k, "region_id"] <- paste(sort(region_info[idx, "region_id"]), collapse = ";")
  }

  # compute LD for weight variants on each chromosome
  chrs <- sort(unique(weight_info$chrom))
  for (b in chrs) {
    loginfo("Computing LD for weight variants on chr%s", b)
    weightinfo <- weight_info[weight_info$chrom == b, ]
    if (nrow(weightinfo) > 0) {
      weight_region_ids <- names(sort(-table(weightinfo$region_id)))
      weight_LD_list <- mclapply(weight_region_ids, function(x){
        # load the R_snp and SNP info for the region
        # and extract LD for the weight variants
        curr_region_LD_list <- list()
        curr_region_ids <- unlist(strsplit(x, ";"))
        curr_region_idx <- match(curr_region_ids, LD_info$region_id)
        LD_matrix_files <- LD_info$LD_matrix[curr_region_idx]
        if (length(LD_matrix_files) > 1) {
          R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader = LD_loader)
          R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
        } else {
          R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader = LD_loader)
        }

        snpinfo <- do.call(rbind, snp_info[curr_region_ids])
        rownames(R_snp) <- snpinfo$id
        colnames(R_snp) <- snpinfo$id

        weight_ids <- weightinfo[weightinfo$region_id == x, "weight_id"]

        for (weight_id in weight_ids) {
          snpnames <- rownames(weights[[weight_id]]$wgt)
          R_wgt <- R_snp[snpnames, snpnames, drop=F]
          curr_region_LD_list[[weight_id]] <- R_wgt
        }
        curr_region_LD_list
      })
      if (length(weight_LD_list) != length(weight_region_ids)) {
        stop("Not all cores returned results. Try rerun with bigger memory or fewer cores")
      }
      weight_LD_list <- unlist(weight_LD_list, recursive = FALSE)
      for(weight_id in names(weight_LD_list)){
        weights[[weight_id]][["R_wgt"]] <- weight_LD_list[[weight_id]]
      }
    }
  }
  return(weights)
}

