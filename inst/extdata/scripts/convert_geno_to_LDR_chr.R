args = commandArgs(trailingOnly=TRUE)

b <- as.numeric(args[1])
reg_start <- as.numeric(args[2])

print(paste0("chromsome: ", b))
print(paste0("start region: ", reg_start))

####################

setwd("/gpfs/data/xhe-lab/ukb_LDR/matrices_0.1")

library(ctwas)
library(tools)
library(data.table)

# specify LD reference
ldref_dir <- "/gpfs/data/xhe-lab/ukb_LDR/genotype_data_0.1"
ldref_prefix <- "ukb_chr"
ldref_suffix <- ""
ldref_filtype <- ".pgen"
ldref_files <- paste0(ldref_dir, "/", ldref_prefix, 1:22, ldref_suffix,ldref_filtype)

# specify LD regions
ld_regions <- "EUR"
ld_regions_custom = NULL

# specify output
out_filestem <- "ukb_b38_0.1_chr"
out_dir <- "/scratch/wcrouse/LDR_b38_cova"

####################

dir.create(out_dir, showWarnings=F)
if (is.null(ld_regions_custom)) {
  regionfile <- system.file("extdata", "ldetect", paste0(ld_regions, ".b38.bed"), package = "ctwas")
} else {
  regionfile <- ld_regions_custom
}
reg <- read.table(regionfile, header = T, stringsAsFactors = F)
if (is.character(reg$chr)) {
  reg$chr <- readr::parse_number(reg$chr)
}
reg <- reg[order(reg$chr, reg$start), ]
ld_pvarfs <- sapply(ldref_files, prep_pvar, outputdir = out_dir)

#load build 38 positions
snpinfo_hg38_all <- fread("/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants_hg38.bim")
colnames(snpinfo_hg38_all) <- c("chr", "id", "cm", "pos", "alt", "ref")

#print(paste0("chromosome ", b))
logging::loginfo(paste0("chromosome ", b))

ld_pvarf <- ld_pvarfs[b]
snpinfo <- read_pvar(ld_pvarf)
if (unique(snpinfo$chrom) != b) {
  stop("Input genotype file not split by chromosome or not in correct order")
}
regions <- reg[reg$chr == b, ]
ld_pgen <- prep_pgen(pgenf = ldref_files[b], ld_pvarfs[b])

snpinfo_hg38 <- snpinfo_hg38_all[snpinfo_hg38_all$chr==b,]
rm(snpinfo_hg38_all)###for single chromosome only

for (rn in reg_start:nrow(regions)) {
  
  #print(rn)
  logging::loginfo(rn)
  
  p0 <- regions[rn, "start"]
  p1 <- regions[rn, "stop"]
  
  outfile <- paste0(out_dir,"/", out_filestem, b, ".R_snp.", p0, "_", p1, ".RDS")
  outfile_Rvar <- paste0(out_dir,"/", out_filestem, b, ".R_snp.", p0, "_", p1, ".Rvar")
  outfile_temp <- paste0(outfile_Rvar, "-temp")
  
  if (!(file.exists(outfile_Rvar)| file.exists(outfile_temp))){
    
    file.create(outfile_temp)
    
    sidx_hg38 <- which(snpinfo_hg38$pos >= p0 & snpinfo_hg38$pos < p1)
    sid_hg38 <- snpinfo_hg38$id[sidx_hg38]
    
    sidx <- which(snpinfo$id %in% sid_hg38)
    sid <- snpinfo[sidx, "id"]
    
    X.g <- read_pgen(ld_pgen, variantidx = sidx)
    
    ptm <- proc.time()
    R_snp <- Rfast::cora(X.g)
    proc.time() - ptm
    logging::loginfo(paste0("runtime: ", (proc.time() - ptm)[3]))
    
    rownames(R_snp) <- sid$id
    colnames(R_snp) <- sid$id
    R_snp <- R_snp[sid_hg38, sid_hg38]
    rownames(R_snp) <- NULL
    colnames(R_snp) <- NULL
    saveRDS(R_snp, file=outfile)
    
    R_snp_variances <- Rfast::colVars(X.g)
    names(R_snp_variances) <- sid$id
    R_snp_variances <- R_snp_variances[sid_hg38]
    names(R_snp_variances) <- NULL

    R_snp_AFs <- colSums(X.g)/(2*nrow(X.g))
    names(R_snp_AFs) <- sid$id
    R_snp_AFs <- R_snp_AFs[sid_hg38]
    names(R_snp_AFs) <- NULL

    R_var <- cbind(snpinfo_hg38[sidx_hg38,c("chr","id","pos","alt","ref")], R_snp_variances, R_snp_AFs)
    colnames(R_var) <- c("chrom","id","pos","alt","ref","variance","allele_freq")
    
    write.table(R_var, file=outfile_Rvar, sep="\t", col.names=T, row.names=F, quote=F)
    
    file.remove(outfile_temp)
  }
}

