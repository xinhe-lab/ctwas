
#' Preprocess PredictDB weights and harmonize with LD reference
#' (adapted from preharmonize_wgt_ld)
#'
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.
#'
#' @param LD_R_dir a string, pointing to a directory containing all LD matrix files and variant information. Expects .RDS files which contain LD correlation matrices for a region/block.
#' For each RDS file, a file with same base name but ended with .Rvar needs to be present in the same folder. the .Rvar file has 5 required columns: "chrom", "id", "pos", "alt", "ref".
#' If using PredictDB format weights and \code{scale_by_ld_variance=T}, a 6th column is also required: "variance", which is the variance of the each SNP.
#' The order of rows needs to match the order of rows in .RDS file.
#'
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name.
#'
#' @param drop_strand_ambig TRUE/FALSE, if TRUE remove strand ambiguous variants (A/T, G/C).
#'
#' @importFrom logging addHandler loginfo
#'
#' @export
#'
preprocess_wgt_ld <- function (weight,
                               ld_R_dir,
                               outputdir = getwd(),
                               outname,
                               drop_strand_ambig = TRUE){

  dir.create(outputdir, showWarnings = F)

  # read the PredictDB weights
  sqlite <- RSQLite::dbDriver("SQLite")
  db = RSQLite::dbConnect(sqlite, weight)
  query <- function(...) RSQLite::dbGetQuery(db, ...)
  weight_table <- query("select * from weights")
  extra_table <- query("select * from extra")
  RSQLite::dbDisconnect(db)

  gnames <- unique(weight_table$gene)
  loginfo("Number of genes with weights provided: %s", length(gnames))

  # load LD information and subset to variants in weights
  ld_Rfs <- write_ld_Rf(ld_R_dir, outname = outname, outputdir = outputdir)

  ld_Rinfo <- lapply(ld_Rfs, data.table::fread, header=T)
  ld_Rinfo <- as.data.frame(do.call(rbind, ld_Rinfo))

  ld_snpinfo <- lapply(ld_Rfs, ctwas:::read_ld_Rvar)
  ld_snpinfo <- as.data.frame(do.call(rbind, ld_snpinfo))

  # remove variants in weight table, but not in LD reference
  loginfo("Number of variants in weights: %s", length(unique(weight_table$rsid)))
  loginfo("Remove %s variants in weights but not in LD reference", length(setdiff(weight_table$rsid, ld_snpinfo$id)))
  weight_table <- weight_table[weight_table$rsid %in% ld_snpinfo$id, ]

  # remove genes with no variants in LD reference
  loginfo("Remove %s genes with no variants in LD reference", length(setdiff(gnames, weight_table$gene)))
  gnames <- unique(weight_table$gene)
  loginfo("Number of genes left after removing variants not in LD reference: %s", length(gnames))

  # subset to variants in weight table
  ld_snpinfo <- ld_snpinfo[ld_snpinfo$id %in% weight_table$rsid,]

  weight_table_harmonized <- list()

  loginfo("Processing weights for %s genes to match LD reference", length(gnames))

  for (i in 1:length(gnames)){

    if (i %% 1000 == 0){
      loginfo("Current gene: %s", i)
    }

    gname <- gnames[i]

    wgt <- weight_table[weight_table$gene==gname,]
    wgt.matrix <- as.matrix(wgt[, "weight", drop = F])

    rsid_varID <- wgt[,c("rsid", "varID")]

    rownames(wgt.matrix) <- wgt$rsid
    chrpos <- do.call(rbind, strsplit(wgt$varID, "_"))

    snps <- data.frame(gsub("chr", "", chrpos[, 1]), wgt$rsid,
                       "0", chrpos[, 2], wgt$eff_allele, wgt$ref_allele,
                       stringsAsFactors = F)
    colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")
    snps$chrom <- as.integer(snps$chrom)
    snps$pos <- as.integer(snps$pos)

    chrom <- unique(snps$chrom)
    ld_Rinfo_chrom <- ld_Rinfo[ld_Rinfo$chrom==chrom,]

    w <- harmonize_wgt_ld(wgt.matrix, snps, ld_snpinfo, drop_strand_ambig)

    wgt.matrix <- w[["wgt"]]
    snps <- w[["snps"]]

    wgt.matrix <- wgt.matrix[abs(wgt.matrix[, "weight"]) > 0, , drop = F]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix),, drop = F]

    snpnames <- intersect(rownames(wgt.matrix), ld_snpinfo$id)

    wgt.idx <- match(snpnames, rownames(wgt.matrix))
    wgt <- wgt.matrix[wgt.idx, "weight", drop = F]

    snps.idx <- match(snpnames, snps$id)
    snps <- snps[snps.idx,]

    if (length(snpnames)>0){
      weight_table_current <- data.frame(gene=gname,
                                         rsid=snps$id,
                                         varID=rsid_varID$varID[match(snps$id, rsid_varID$rsid)],
                                         ref_allele=snps$ref,
                                         eff_allele=snps$alt,
                                         weight=wgt[,"weight"])

      weight_table_harmonized[[gname]] <- weight_table_current
    }
  }
  weight_table_harmonized <- do.call(rbind, weight_table_harmonized)

  loginfo("Number of genes with weights harmonized with LD reference: %s", length(unique(weight_table_harmonized$gene)))

  extra_table <- extra_table[extra_table$gene %in% weight_table_harmonized$gene,]

  if (file.exists(paste0(outputdir, outname, ".db"))){
    invisible(file.remove(paste0(outputdir, outname, ".db")))
  }

  db <- RSQLite::dbConnect(sqlite, file.path(outputdir, paste0(outname, ".db")))
  RSQLite::dbWriteTable(db, "extra", extra_table)
  RSQLite::dbWriteTable(db, "weights", weight_table_harmonized)
  RSQLite::dbDisconnect(db)

  invisible(file.remove(file.path(outputdir, paste0(outname, "_ld_R_chr", 1:22, ".txt"))))

  loginfo("Save processed weights to: %s", file.path(outputdir, paste0(outname, ".db")))

  return(list(weight_table = weight_table_harmonized, extra_table = extra_table))
}

