## Libraries
# library(ggplot2)
library(logging)
library(foreach)

library(ctwas)
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
outputdir <- paste0("/project2/xinhe/shared_data/multigroup_ctwas/test_data/output/", trait, ".", tissue, "/")
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
processed_z_snp_file <- file.path(outputdir, paste0(outname, ".preprocessed.z_snp.Rd"))
if (file.exists(processed_z_snp_file)){
  cat(sprintf("Load preprocessed z_snp: %s \n", processed_z_snp_file))
  load(processed_z_snp_file)
}else{
  load(z_snp_outfile)
  runtime <- system.time({
    z_snp <- preprocess_z_snp(z_snp,
                              region_info,
                              gwas_n,
                              drop_multiallelic = TRUE,
                              drop_strand_ambig = TRUE)
    save(z_snp, file = file.path(outputdir, paste0(outname, ".preprocessed.z_snp.Rd")))
  })
  cat(sprintf("Preprocessing GWAS z-scores took %0.2f minutes\n",runtime["elapsed"]/60))
}

## Preprocess weights
processed_weight_file <- file.path(outputdir, paste0(outname, ".preprocessed.weights.Rd"))
if (file.exists(processed_weight_file)){
  cat(sprintf("Load preprocessed weight: %s\n", processed_weight_file))
  load(processed_weight_file)
}else{
  runtime <- system.time({
    weights <- preprocess_weights(weight_file,
                                  region_info,
                                  z_snp,
                                  drop_strand_ambig = TRUE,
                                  filter_protein_coding_genes = TRUE,
                                  scale_by_ld_variance = TRUE)
    save(weights, file = processed_weight_file)
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
    z_gene <- compute_gene_z(z_snp, weights, logfile = file.path(outputdir, paste0(outname, ".compute_gene_z.log")))
    save(z_gene, file = gene_z_file)
  })
  cat(sprintf("Imputing gene z-scores took %0.2f minutes\n",runtime["elapsed"]/60))
}

# combine z-scores of SNPs and genes
old_outputdir <- paste0("/project2/xinhe/shared_data/multigroup_ctwas/test_data/output/", trait, "_", tissue, "/")
processed_z_snp_file <- file.path(outputdir, paste0(outname, ".preprocessed.z_snp.Rd"))
processed_weight_file <- file.path(outputdir, paste0(outname, ".preprocessed.weights.Rd"))
gene_z_file <- file.path(old_outputdir, paste0(outname, ".gene_z.Rd"))
load(processed_z_snp_file)
load(processed_weight_file)
load(gene_z_file)

# get regionlist
regionlist_thin_file <- file.path(outputdir, paste0(outname, ".regionlist.thin", thin, ".RDS"))
if (file.exists(regionlist_thin_file)) {
  res <- readRDS(regionlist_thin_file)
} else{
  loginfo("Get regionlist with thin = %.2f", thin)
  res <- get_regionlist(region_info,
                        z_snp,
                        z_gene,
                        weights,
                        maxSNP = max_snp_region,
                        trim_by = "random",
                        thin = thin,
                        minvar = 2,
                        mingene = 0,
                        adjust_boundary_genes = TRUE)
  saveRDS(res, regionlist_thin_file)
}
regionlist <- res$regionlist
weights <- res$weights
boundary_genes <- res$boundary_genes
rm(res)

# old_regionlist_res <- readRDS(file.path(old_outputdir, paste0(outname, ".regionlist.thin", thin, ".RDS")))
# identical(regionlist[[100]]$gid, old_regionlist_res$regionlist[[100]]$gid)
# identical(regionlist[[100]]$sid, old_regionlist_res$regionlist[[100]]$sid)

##### ctwas_rss parameter estimation #####
cat("##### Estimate parameters ##### \n")
## est_param
param_file <- file.path(outputdir, paste0(outname, ".param.RDS"))
if (file.exists(param_file)) {
  param <- readRDS(param_file)
} else{
  runtime <- system.time({
    param <- est_param(z_snp,
                       z_gene,
                       regionlist,
                       thin = thin,
                       niter1 = 3,
                       niter2 = 30,
                       group_prior_var_structure = "independent",
                       logfile = file.path(outputdir, paste0(outname, ".est_param.log")),
                       ncore = 6)
  })
  cat(sprintf("Parameter estimation took %0.2f minutes\n",runtime["elapsed"]/60))
  saveRDS(param, param_file)
}
group_prior <- param$group_prior
group_prior_var <- param$group_prior_var

old_param <- readRDS(file.path(old_outputdir, paste0(outname, ".param.RDS")))

# Assessing parameter estimates
ctwas_parameters <- summarize_param(param, gwas_n)
saveRDS(ctwas_parameters, paste0(outputdir, "/", outname, ".ctwas_parameters.RDS"))

##### Screen regions #####
load("/project/xinhe/shengqian/cTWAS_analysis/data/make_test_data/make_test_data_mergeoff.s2.susieIrssres.Rd")
group_prior <- group_prior_rec[,ncol(group_prior_rec)]
group_prior_var <- group_prior_var_rec[,ncol(group_prior_var_rec)]
group_prior["SNP"] <- group_prior["SNP"] * thin # adjust parameter to account for thin argument

runtime <- system.time({
  screened_region_tags <- screen_regions(z_snp,
                                         z_gene,
                                         regionlist,
                                         region_info,
                                         weights,
                                         thin = thin,
                                         L = 5,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         max_snp_region = max_snp_region,
                                         ncore = ncore,
                                         logfile = file.path(outputdir, paste0(outname, ".screen_regions.L5.log")),
                                         verbose = TRUE)
})
saveRDS(screened_region_tags, file.path(outputdir, paste0(outname, ".screen_regions.L5.res.RDS")))
cat(sprintf("Screen regions took %0.2f minutes\n",runtime["elapsed"]/60))

# Expand screened regionlist with all SNPs in the regions
screened_regionlist <- regionlist[screened_region_tags]
if (thin < 1){
  loginfo("Update regionlist with full SNPs for screened regions")
  screened_regionlist <- expand_regionlist(screened_regionlist,
                                           region_info,
                                           z_snp,
                                           trim_by = "z",
                                           maxSNP = max_snp_region)
}

##### finemapping #####
res <- readRDS(file.path(old_outputdir, paste0(outname, ".screen_regions.L5.max_iter100.res.RDS")))
old_screened_region_tags <- res$screened_region_tags
old_screened_regionlist <- res$screened_regionlist
rm(res)

screened_regionlist <- old_screened_regionlist
# if (setequal(screened_region_tags, old_screened_region_tags)){
#   cat("screened_regions PASS")
# }else {
#   cat("screened_regions FAIL")
# }

# Finemap a single region
region_tag <- "16:71020125-72901251"
runtime <- system.time({
  finemap_region_res <- finemap_region(z_snp,
                                       z_gene,
                                       region_tag = region_tag,
                                       region_info = region_info,
                                       weights = weights,
                                       L = 5,
                                       group_prior = group_prior,
                                       group_prior_var = group_prior_var,
                                       force_compute_cor = TRUE,
                                       save_cor = TRUE,
                                       cor_dir = file.path(outputdir, "cor_matrix"),
                                       verbose = TRUE)
})
cat(sprintf("Finemapping region took %0.2f seconds \n",runtime["elapsed"]))

# old_finemap_full_res <- as.data.frame(data.table::fread("/project/xinhe/shengqian/cTWAS_analysis/data/make_test_data/make_test_data_mergeoff.susieIrss.txt", header = T))
# old_finemap_full_res$id[old_finemap_full_res$type == "gene"] <-
#   paste0(old_finemap_full_res$id[old_finemap_full_res$type == "gene"], "|", "mashr_Liver_nolnc")
#
# ids <- intersect(old_finemap_full_res$id, finemap_region_res$id)
# old_finemap_region_res <- old_finemap_full_res[match(ids, old_finemap_full_res$id), ]
# new_finemap_region_res <- finemap_region_res[match(ids, finemap_region_res$id), ]
#
# plot(old_finemap_region_res$susie_pip, new_finemap_region_res$susie_pip)
#
# all.equal(old_finemap_region_res$susie_pip, new_finemap_region_res$susie_pip)
# all.equal(old_finemap_region_res$mu2, new_finemap_region_res$mu2)
# all.equal(old_finemap_region_res$cs_index, new_finemap_region_res$cs_index)

## Finemapping screened regions
cat("##### Finemapping screened regions ##### \n")
region_tag <- "16:71020125-72901251"
finemap_region_res <- finemap_region(z_snp,
                                     z_gene,
                                     region_tag = region_tag,
                                     region_info = region_info,
                                     weights = weights,
                                     L = 5,
                                     group_prior = group_prior,
                                     group_prior_var = group_prior_var,
                                     force_compute_cor = TRUE,
                                     save_cor = TRUE,
                                     cor_dir = file.path(outputdir, "cor_matrix"),
                                     verbose = TRUE)

runtime <- system.time({
  finemap_res <- finemap_regions(z_snp,
                                 z_gene,
                                 screened_regionlist,
                                 region_info,
                                 weights,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 L = 5,
                                 force_compute_cor = TRUE,
                                 save_cor = TRUE,
                                 cor_dir = paste0(outputdir, "/cor_matrix/"),
                                 ncore = ncore,
                                 logfile = file.path(outputdir, paste0(outname, ".finemapping_regions.log")),
                                 verbose = TRUE)
})
saveRDS(finemap_res, file.path(outputdir, paste0(outname, ".finemap_regions.res.RDS")))
cat(sprintf("Finemapping took %0.2f minutes\n",runtime["elapsed"]/60))

# Print sessionInfo
sessionInfo()
