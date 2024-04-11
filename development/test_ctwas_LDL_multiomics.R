## Libraries
library(ctwas)
library(ggplot2)
library(tools)
library(logging)

##### Settings #####
trait <- "LDL"
tissue <- "Multiomics"
gwas_file <- "/project2/xinhe/shared_data/multigroup_ctwas/test_data/ukb-d-30780_irnt.vcf.gz"
gwas_n <- 343621
genome_version <- "b38"
ld_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1/"
weight_files <- c("/project/xinhe/shengqian/cTWAS_analysis/data/Liver_Expression.db",
                 "/project/xinhe/shengqian/cTWAS_analysis/data/Lung_Expression.db",
                 "/project/xinhe/shengqian/cTWAS_analysis/data/Liver_Splicing.db",
                 "/project/xinhe/shengqian/cTWAS_analysis/data/Lung_Splicing.db")
thin <- 0.1
max_snp_region <- 20000
ncore <- 6
outputdir <- "/project2/xinhe/shared_data/multigroup_ctwas/test_data/output/LDL_Multiomics_Kevin/"
outname <- paste0(trait, ".", tissue)
dir.create(outputdir, showWarnings=F, recursive=T)

##### LD region info #####
region_info_file <- file.path(outputdir, "region_info.txt")
if (file.exists(region_info_file)){
  region_info <- read.table(region_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
} else {
  region_file <- system.file("extdata/ldetect", paste0("EUR.", genome_version, ".bed"), package = "ctwas")
  region_info <- read.table(region_file, header = TRUE)
  colnames(region_info)[1:3] <- c("chrom", "start", "stop")
  region_info$chrom <- as.numeric(gsub("chr", "", region_info$chrom))
  region_info$region_id <- paste0(region_info$chr, ":", region_info$start, "-", region_info$stop)
  filestem <- paste0("ukb_", genome_version, "_0.1")
  ld_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_info$chrom, region_info$start, region_info$stop)
  region_info$LD_matrix <- file.path(ld_R_dir, paste0(ld_filestem, ".RDS"))
  region_info$SNP_info <- file.path(ld_R_dir, paste0(ld_filestem, ".Rvar"))
  stopifnot(all(file.exists(region_info$LD_matrix)))
  stopifnot(all(file.exists(region_info$SNP_info)))
  write.table(region_info, file = file.path(outputdir, "region_info.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
}

##### Preprocess GWAS z-scores #####
cat("##### Preprocess z-scores ##### \n")
gwas_name <- "ukb-d-30780_irnt"
z_snp_outfile <- file.path(outputdir, paste0(gwas_name, ".z_snp.Rd"))
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

  # save the formatted z-scores and GWAS sample size
  save(z_snp, gwas_n, file=z_snp_outfile)
}

## Preprocess z-scores
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
  })
  save(z_snp, file = file.path(outputdir, paste0(outname, ".preprocessed.z_snp.Rd")))
  cat(sprintf("Preprocessing GWAS z-scores took %0.2f minutes\n",runtime["elapsed"]/60))
}

##### Preprocess weights #####
processed_weight_file <- file.path(outputdir, paste0(outname, ".preprocessed.weights.Rd"))
if (file.exists(processed_weight_file)){
  cat(sprintf("Load preprocessed weight: %s\n", processed_weight_file))
  load(processed_weight_file)
}else{
  runtime <- system.time({
    weights_liver_expression <- preprocess_weights(weight_file = weight_files[1],
                                                   region_info,
                                                   z_snp,
                                                   type = "expression",
                                                   context = "liver",
                                                   weight_format = "PredictDB",
                                                   ncore = ncore,
                                                   drop_strand_ambig = TRUE,
                                                   scale_by_ld_variance = TRUE,
                                                   load_predictdb_LD = TRUE,
                                                   filter_protein_coding_genes = TRUE)

    weights_lung_expression <- preprocess_weights(weight_file = weight_files[2],
                                                  region_info,
                                                  z_snp,
                                                  type = "expression",
                                                  context = "lung",
                                                  weight_format = "PredictDB",
                                                  ncore = ncore,
                                                  drop_strand_ambig = TRUE,
                                                  scale_by_ld_variance = TRUE,
                                                  load_predictdb_LD = TRUE,
                                                  filter_protein_coding_genes = TRUE)

    weights_liver_splicing <- preprocess_weights(weight_file = weight_files[3],
                                                 region_info,
                                                 z_snp,
                                                 type = "splicing",
                                                 context = "liver",
                                                 weight_format = "PredictDB",
                                                 ncore = ncore,
                                                 drop_strand_ambig = TRUE,
                                                 scale_by_ld_variance = TRUE,
                                                 load_predictdb_LD = TRUE,
                                                 filter_protein_coding_genes = FALSE)

    weights_lung_splicing <- preprocess_weights(weight_file = weight_files[4],
                                                region_info,
                                                z_snp,
                                                type = "splicing",
                                                context = "lung",
                                                weight_format = "PredictDB",
                                                ncore = ncore,
                                                drop_strand_ambig = TRUE,
                                                scale_by_ld_variance = TRUE,
                                                load_predictdb_LD = TRUE,
                                                filter_protein_coding_genes = FALSE)

    weights <- c(weights_liver_expression, weights_lung_expression, weights_liver_splicing, weights_lung_splicing)
  })
  save(weights, file = processed_weight_file)
  cat(sprintf("Preprocessing weights took %0.2f minutes\n",runtime["elapsed"]/60))
}

##### Impute gene z-scores #####
cat("##### Imputing gene z-scores ##### \n")
gene_z_file <- file.path(outputdir, paste0(outname, ".gene_z.Rd"))
if( file.exists(gene_z_file) ){
  cat(sprintf("Load gene z-scores from %s \n", gene_z_file))
  load(gene_z_file)
}else{
  runtime <- system.time({
    z_gene <- compute_gene_z(z_snp, weights, ncore=ncore,
                             logfile = file.path(outputdir, paste0(outname, ".compute_gene_z.log")))
  })
  save(z_gene, file = gene_z_file)
  cat(sprintf("Imputing gene z-scores took %0.2f minutes\n",runtime["elapsed"]/60))
}

##### Get regionlist #####
regionlist_thin_file <- file.path(outputdir, paste0(outname, ".regionlist.thin", thin, ".RDS"))
if (file.exists(regionlist_thin_file)) {
  res <- readRDS(regionlist_thin_file)
} else{
  runtime <- system.time({
    loginfo("Get regionlist with thin = %.2f", thin)
    res <- get_regionlist(region_info,
                          z_snp,
                          z_gene,
                          weights,
                          maxSNP = max_snp_region,
                          trim_by = "random",
                          thin = thin,
                          minvar = 2,
                          adjust_boundary_genes = TRUE,
                          ncore = ncore)
  })
  saveRDS(res, regionlist_thin_file)
  cat(sprintf("Get regionlist took %0.2f minutes\n",runtime["elapsed"]/60))
}
regionlist <- res$regionlist
boundary_genes <- res$boundary_genes
rm(res)

##### Estimate parameters #####
cat("##### Estimating parameters ##### \n")
## est_param
param_file <- file.path(outputdir, paste0(outname, ".param.RDS"))
if (file.exists(param_file)) {
  param <- readRDS(param_file)
} else{
  runtime <- system.time({
    param <- est_param(regionlist,
                       niter1 = 3,
                       niter2 = 30,
                       group_prior_var_structure = "independent",
                       logfile = file.path(outputdir, paste0(outname, ".est_param.log")),
                       ncore = ncore)
  })
  saveRDS(param, param_file)
  cat(sprintf("Parameter estimation took %0.2f minutes\n",runtime["elapsed"]/60))
}
group_prior <- param$group_prior
group_prior_var <- param$group_prior_var


##### Assess parameter estimates #####
ctwas_parameters <- summarize_param(param, gwas_n)
saveRDS(ctwas_parameters, paste0(outputdir, "/", outname, ".ctwas_parameters.RDS"))

##### Screen regions #####
screen_regions_file <- file.path(outputdir, paste0(outname, ".screened_regionlist.RDS"))
if (file.exists(screen_regions_file)) {
  screened_regionlist <- readRDS(screen_regions_file)
} else{
  runtime <- system.time({
    screened_region_ids <- screen_regions(regionlist,
                                           region_info,
                                           weights,
                                           L = 5,
                                           group_prior = group_prior,
                                           group_prior_var = group_prior_var,
                                           max_snp_region = max_snp_region,
                                           ncore = ncore,
                                           verbose = TRUE,
                                           logfile = file.path(outputdir, paste0(outname, ".screen_regions.log")))
  })
  cat(sprintf("Screen regions took %0.2f minutes\n",runtime["elapsed"]/60))
  loginfo("%d regions left after screening regions", length(screened_region_ids))

  # Expand screened regionlist with all SNPs in the regions
  screened_regionlist <- regionlist[screened_region_ids]
  if (thin < 1){
    loginfo("Expand regionlist with full SNPs for %d screened regions", length(screened_regionlist))
    screened_regionlist <- expand_regionlist(screened_regionlist,
                                             region_info,
                                             z_snp,
                                             z_gene,
                                             trim_by = "z",
                                             maxSNP = max_snp_region,
                                             ncore = ncore)
  }
  saveRDS(screened_regionlist, screen_regions_file)
}

##### Finemapping #####

## Finemapping screened regions
cat("##### Finemapping screened regions ##### \n")
runtime <- system.time({
  finemap_res <- finemap_regions(screened_regionlist,
                                 region_info,
                                 weights,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 L = 5,
                                 force_compute_cor = TRUE,
                                 save_cor = TRUE,
                                 cor_dir = paste0(outputdir, "/cor_matrix/"),
                                 ncore = ncore,
                                 verbose = TRUE)
})
saveRDS(finemap_res, file.path(outputdir, paste0(outname, ".finemap_regions.res.RDS")))
cat(sprintf("Finemapping took %0.2f minutes\n",runtime["elapsed"]/60))
