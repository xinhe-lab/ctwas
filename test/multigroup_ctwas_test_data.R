## Libraries
# library(ctwas)
# library(ggplot2)
devtools::load_all("/home/kaixuan/projects/cTWAS_package/multigroup_test/ctwas/.")

## Settings
trait <- "LDL"
tissue <- "Liver"
datadir <- "/project2/xinhe/shared_data/multigroup_ctwas/test_data/"
weight_file <- "/project2/xinhe/shared_data/multigroup_ctwas/test_data/mashr_Liver_nolnc.db"
gwas_file <- "/project2/xinhe/shared_data/multigroup_ctwas/test_data/ukb-d-30780_irnt.vcf.gz"
gwas_n <- 343621
ncore <- 4
thin <- 0.1
max_snp_region <- 20000
outputdir <- paste0("/project2/xinhe/shared_data/multigroup_ctwas/test_data/output/", trait, "_", tissue, "/")
outname <- paste0(trait, ".", tissue)
dir.create(outputdir, showWarnings=F, recursive=T)

## LD reference and region info
region_info_file <- file.path(outputdir, "region_info.txt")
if (file.exists(region_info_file)){
  region_info <- read.table(region_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
} else {
  region_file <- system.file("extdata/ldetect", "EUR.b38.bed", package = "ctwas")
  region_info <- read.table(region_file, header = TRUE)
  colnames(region_info)[1:3] <- c("chrom", "start", "stop")
  region_info$chrom <- as.numeric(gsub("chr", "", region_info$chrom))
  region_info$region_tag <- paste0(region_info$chr, ":", region_info$start, "-", region_info$stop)

  ld_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1"
  filestem <- "ukb_b38_0.1"
  ld_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_info$chrom, region_info$start, region_info$stop)
  region_info$LD_matrix <- file.path(ld_R_dir, paste0(ld_filestem, ".RDS"))
  region_info$SNP_info <- file.path(ld_R_dir, paste0(ld_filestem, ".Rvar"))
  write.table(region_info, file = file.path(outputdir, "region_info.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
}

## Prepare GWAS z-scores
gwas_name <- "ukb-d-30780_irnt"
z_snp_outfile <- file.path(outputdir, paste0(gwas_name, ".z_snp.Rd"))
print(z_snp_outfile)

if (!file.exists(z_snp_outfile)){
  # read the data using the VariantAnnotation package
  z_snp <- VariantAnnotation::readVcf("/project2/xinhe/shared_data/multigroup_ctwas/test_data/ukb-d-30780_irnt.vcf.gz")
  z_snp <- as.data.frame(gwasvcf::vcf_to_tibble(z_snp))

  # compute the z-scores
  z_snp$Z <- z_snp$ES/z_snp$SE

  # collect sample size (most frequent sample size for all variants)
  gwas_n <- as.numeric(names(sort(table(z_snp$SS),decreasing=TRUE)[1]))
  cat("gwas_n=", gwas_n, "\n")

  # subset the columns and format the column names
  z_snp <- z_snp[,c("rsid", "ALT", "REF", "Z")]
  colnames(z_snp) <- c("id", "A1", "A2", "z")

  # # drop multiallelic variants (id not unique)
  # z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),]

  # save the formatted z-scores and GWAS sample size
  save(z_snp, gwas_n, file=z_snp_outfile)
}

## Harmonize z-scores with LD reference
processed_z_snp_file <- file.path(outputdir, paste0(outname, ".processed.z_snp.Rd"))
if (file.exists(processed_z_snp_file)){
  cat(sprintf("Load preprocessed z_snp: %s \n", processed_z_snp_file))
  load(processed_z_snp_file)
}else{
  load(z_snp_outfile)
  runtime <- system.time({
    res <- process_z(z_snp,
                     region_info,
                     gwas_n,
                     drop_multiallelic = TRUE,
                     drop_strand_ambig = TRUE,
                     detect_ld_mismatch = FALSE,
                     ncore = ncore)
    z_snp <- res$z_snp
    ld_mismatch_res <- res$ld_mismatch_res
    save(z_snp, ld_mismatch_res, file = file.path(outputdir, paste0(outname, ".processed.z_snp.Rd")))
  })
  cat(sprintf("Preprocessing GWAS z-scores took %0.2f minutes\n",runtime["elapsed"]/60))

  # runtime <- system.time({
  #   res <- process_z(z_snp,
  #                    region_info,
  #                    gwas_n,
  #                    drop_multiallelic = TRUE,
  #                    drop_strand_ambig = TRUE,
  #                    detect_ld_mismatch = TRUE,
  #                    ncore = ncore)
  #   z_snp <- res$z_snp
  #   ld_mismatch_res <- res$ld_mismatch_res
  #   save(z_snp, ld_mismatch_res, file = file.path(outputdir, paste0(outname, ".processed.z_snp.ld_mismatch.Rd")))
  # })
  # cat(sprintf("Preprocessing GWAS z-scores took %0.2f minutes\n",runtime["elapsed"]/60))
}

# # compare with earlier version
# z_snp_new <- z_snp
# z_snp_old <- get(load("/project/xinhe/shengqian/cTWAS_analysis/data/make_test_data/make_test_data_z_snp.Rd"))
# z_snp_old <- z_snp_old[, c("id", "A1", "A2", "z")]
# if (all.equal(z_snp_new, z_snp_old)){
#   cat("z_snp check: PASS\n")
# }else{
#   cat("z_snp check: FAIL\n")
# }
# z_snp <- z_snp_new
# rm(z_snp_new, z_snp_old)

## Preprocess weights
processed_weight_file <- file.path(outputdir, paste0(outname, ".processed.weights.Rd"))
if (file.exists(processed_weight_file)){
  cat(sprintf("Load preprocessed weight: %s\n", processed_weight_file))
  load(processed_weight_file)
}else{
  runtime <- system.time({
    res <- process_weights(weight_file,
                           region_info,
                           z_snp,
                           drop_strand_ambig = TRUE,
                           filter_protein_coding_genes = TRUE,
                           scale_by_ld_variance = TRUE)
    weight_list <- res$weight_list
    weight_info <- res$weight_info
    save(weight_list, weight_info, file = processed_weight_file)
  })
  cat(sprintf("Preprocessing weights took %0.2f minutes\n",runtime["elapsed"]/60))
}

## Imputing gene z-scores
cat("##### Imputing gene z-scores ##### \n")
gene_z_file <- file.path(outputdir, paste0(outname, ".gene_z.Rd"))
if( file.exists(gene_z_file) ){
  cat(sprintf("Load gene z-scores from %s \n", gene_z_file))
  load(gene_z_file)
}else{
  runtime <- system.time({
    res <- compute_gene_z(z_snp, region_info, weight_list, weight_info,
                          logfile = file.path(outputdir, paste0(outname, ".compute_gene_z.log")),
                          ncore=ncore)
    z_gene <- res$z_gene
    gene_info <- res$gene_info
    save(z_gene, gene_info, file = gene_z_file)
  })
  cat(sprintf("Imputing gene z-scores took %0.2f minutes\n",runtime["elapsed"]/60))
}

# # compare with earlier version
# load(gene_z_file)
# z_gene_new <- z_gene
# gene_info_new <- gene_info
#
# load("/project/xinhe/shengqian/cTWAS_analysis/data/make_test_data/make_test_data_z_gene.Rd")
# z_gene_old <- z_gene
# z_gene_old$id <- paste0(z_gene_old$id, "|", "mashr_Liver_nolnc")
# nrow(z_gene_old) == nrow(z_gene_new)
# idx <- match(z_gene_old$id, z_gene_new$id)
# z_gene_new <- z_gene_new[idx, ]
#
# if (all.equal(z_gene_new$z, z_gene_old$z)){
#   cat("z_gene check: PASS\n")
# }else{
#   cat("z_gene check: FAIL\n")
# }
#
# gene_info_old <- as.data.frame(data.table::fread("/project/xinhe/shengqian/cTWAS_analysis/data/make_test_data/make_test_data_chr21.exprvar"))
# gene_info_old$id <- paste0(gene_info_old$gene_name, "|", gene_info_old$weight_name)
# nrow(gene_info_old) == nrow(gene_info_new)
# idx <- match(gene_info_old$id, gene_info_new$id)
# gene_info_new <- gene_info_new[idx, ]
#
# if (all.equal(gene_info_old$gene_name, gene_info_new$gene_name)){
#   cat("gene_info check: PASS\n")
# }else{
#   cat("gene_info check: FAIL\n")
# }

# combine z-scores of SNPs and genes
zdf <- combine_z(z_snp, z_gene)
processed_zdf_file <- file.path(outputdir, paste0(outname, ".processed.zdf.Rd"))
save(zdf, file = processed_zdf_file)

# combine z-scores of SNPs and genes
processed_z_snp_file <- file.path(outputdir, paste0(outname, ".processed.z_snp.Rd"))
gene_z_file <- file.path(outputdir, paste0(outname, ".gene_z.Rd"))
processed_weight_file <- file.path(outputdir, paste0(outname, ".processed.weights.Rd"))
processed_zdf_file <- file.path(outputdir, paste0(outname, ".processed.zdf.Rd"))

load(processed_z_snp_file)
load(gene_z_file)
load(processed_weight_file)
load(processed_zdf_file)

# get regionlist
regionlist_thin_file <- file.path(outputdir, paste0(outname, ".regionlist.thin", thin, ".RDS"))

if (file.exists(regionlist_thin_file)) {
  res <- readRDS(regionlist_thin_file)
} else{
  loginfo("Get regionlist with thin = %.2f", thin)
  res <- get_regionlist(region_info,
                        gene_info,
                        weight_list = weight_list,
                        select = zdf$id,
                        thin = thin,
                        maxSNP = max_snp_region,
                        minvar = 2,
                        adjust_boundary = TRUE)
  saveRDS(res, regionlist_thin_file)
}
regionlist <- res$regionlist
weight_list <- res$weight_list
boundary_genes <- res$boundary_genes
rm(res)

##### ctwas_rss parameter estimation #####
cat("##### Estimate parameters ##### \n")
## est_param
param_file <- file.path(outputdir, paste0(outname, ".param.RDS"))
if (file.exists(param_file)) {
  param <- readRDS(param_file)
} else{
  runtime <- system.time({
    param <- est_param(z_snp = z_snp,
                       z_gene = z_gene,
                       region_info = region_info,
                       gene_info = gene_info,
                       regionlist = regionlist,
                       weight_list = weight_list,
                       thin = thin,
                       max_snp_region = max_snp_region,
                       group_prior_var_structure = "independent",
                       max_iter = 1,
                       logfile = file.path(outputdir, paste0(outname, ".est_param.max_iter1.log")),
                       ncore = ncore)
  })
  saveRDS(param, param_file)
  cat(sprintf("Parameter estimation took %0.2f minutes\n",runtime["elapsed"]/60))
}
summary(param)
group_prior <- param$group_prior
group_prior_var <- param$group_prior_var

identical(regionlist, param$regionlist)
identical(weight_list, param$weight_list)
identical(boundary_genes, param$boundary_genes)
regionlist <- param$regionlist
weight_list <- param$weight_list
boundary_genes <- param$boundary_genes

# Assessing parameter estimates
ctwas_parameters <- summarize_param(param, gwas_n, thin = thin)
saveRDS(ctwas_parameters, paste0(outputdir, "/", outname, ".ctwas_parameters.RDS"))

##### Screen regions #####
file.remove(file.path(outputdir, paste0(outname, ".screen_regions.L5.max_iter100.log")))
runtime <- system.time({
  res <- screen_regions(z_snp,
                        z_gene,
                        region_info = region_info,
                        gene_info = gene_info,
                        weight_list = weight_list,
                        regionlist = regionlist,
                        thin = thin,
                        max_snp_region = max_snp_region,
                        L = 5,
                        group_prior = group_prior,
                        group_prior_var = group_prior_var,
                        max_iter = 100,
                        ncore = ncore,
                        logfile = file.path(outputdir, paste0(outname, ".screen_regions.L5.max_iter100.log")),
                        verbose = TRUE)
})
saveRDS(res, file.path(outputdir, paste0(outname, ".screen_regions.L5.max_iter100.res.RDS")))
cat(sprintf("Screen regions took %0.2f minutes\n",runtime["elapsed"]/60))

screened_regionlist <- res$screened_regionlist
screened_region_tags <- res$screened_region_tags
weak_region_finemap_res <- res$weak_region_finemap_res
rm(res)

# 2024-03-20 10:10:10 INFO::Number of region tags that contain strong gene signals: 60
# Screen regions took 81.52 minutes

##### finemapping #####
cat("##### cTWAS finemapping ##### \n")

# Finemap a single region
region_tag <- "16:71020125-72901251"
region_tag <- "1:10583-1961168"
runtime <- system.time({
  ctwas_res <- finemap_region(z_snp = z_snp,
                              z_gene = z_gene,
                              gene_info = gene_info,
                              region_tag = region_tag,
                              regionlist = regionlist,
                              region_info = region_info,
                              weight_list = weight_list,
                              L = 5,
                              group_prior = group_prior,
                              group_prior_var = group_prior_var,
                              verbose = TRUE)
})
cat(sprintf("Finemapping region took %0.2f seconds \n",runtime["elapsed"]))

cat("Run cTWAS finemapping with L = 5: \n")
# convert SNP prior
group_prior["SNP"] <- group_prior["SNP"] * thin

runtime <- system.time({
  strong_region_finemap_res <- list()
  for (region_tag in screened_region_tags) {
    strong_region_finemap_res[[region_tag]] <- finemap_region(z_snp = z_snp,
                                                              z_gene = z_gene,
                                                              gene_info = gene_info,
                                                              regionlist = screened_regionlist,
                                                              region_tag = region_tag,
                                                              weight_list = weight_list,
                                                              L = 5,
                                                              group_prior = group_prior,
                                                              group_prior_var = group_prior_var,
                                                              save_cor = TRUE,
                                                              cor_dir = paste0(outputdir, "/cor_matrix/"))
  }
  strong_region_finemap_res <- do.call(rbind, strong_region_finemap_res)
})
cat(sprintf("Finemapping took %0.2f minutes\n",runtime["elapsed"]/60))

# Print sessionInfo
sessionInfo()
