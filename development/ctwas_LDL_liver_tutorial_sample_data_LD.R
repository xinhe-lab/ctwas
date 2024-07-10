## Libraries
library(ctwas)
library(logging)
library(data.table)
devtools::load_all("/home/kaixuan/projects/cTWAS_package/single_group/ctwas/")

##### Settings #####
gwas_name <- "ukb-d-30780_irnt"
gwas_file <- "/project2/xinhe/shared_data/multigroup_ctwas/tutorial/LDL_multitissue_tutorial/ukb-d-30780_irnt.vcf.gz"
gwas_n <- 343621
ld_R_dir <- "/project2/mstephens/wcrouse/UKB_LDR_0.1/"
weight_file <- "/project2/xinhe/shared_data/multigroup_ctwas/weights/expression_models/expression_Liver.db"
thin <- 0.1
max_snp_region <- 20000
min_nonSNP_PIP <- 0.5
ncore <- 6
outputdir <- "/project2/xinhe/shared_data/singlegroup_ctwas/tutorial/LDL_liver_tutorial/sample_data/LDL_liver_chr16_example"
dir.create(outputdir, showWarnings=F, recursive=T)
cor_dir <- file.path(outputdir, "/cor_matrix")
outname <- "LDL_example"
example_chrom <- 16

multigroup_outputdir <- "/project2/xinhe/shared_data/multigroup_ctwas/tutorial/LDL_multitissue_tutorial/sample_data/LDL_liver_chr16_noLD_example"

##### LD region info #####
region_info_file <- file.path(outputdir, paste0(outname, ".region_info.RDS"))
snp_info_file <- file.path(outputdir, paste0(outname, ".snp_info.RDS"))
LD_info_file <- file.path(outputdir, paste0(outname, ".LD_info.RDS"))

if (file.exists(region_info_file) && file.exists(snp_info_file) && file.exists(LD_info_file) ){
  cat(sprintf("Load preprocessed region_info: %s \n", region_info_file))
  region_info <- readRDS(region_info_file)
  snp_info <- readRDS(snp_info_file)
  LD_info <- readRDS(LD_info_file)
}else{
  region_file <- system.file("extdata/ldetect", "EUR.b38.ldetect.regions.RDS", package = "ctwas")
  region_info <- readRDS(region_file)

  filestem <- paste0("ukb_b38_0.1")
  ld_filestem <- sprintf("%s_chr%s.R_snp.%s_%s", filestem, region_info$chrom, region_info$start, region_info$stop)
  region_info$LD_matrix <- file.path(ld_R_dir, paste0(ld_filestem, ".RDS"))
  region_info$SNP_info <- file.path(ld_R_dir, paste0(ld_filestem, ".Rvar"))
  res <- preprocess_region_LD_snp_info(region_info,
                                       chrom = example_chrom,
                                       use_LD = TRUE,
                                       ncore = ncore)
  region_info <- res$region_info
  snp_info <- res$snp_info
  LD_info <- res$LD_info

  saveRDS(region_info, region_info_file)
  saveRDS(snp_info, snp_info_file)
  saveRDS(LD_info, LD_info_file)
}

region_info_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".region_info.RDS")))
snp_info_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".snp_info.RDS")))
all.equal(region_info, region_info_multigroup)
all.equal(snp_info, snp_info_multigroup)

##### Preprocess GWAS z-scores #####
z_snp_outfile <- file.path(outputdir, paste0(gwas_name, ".z_snp.RDS"))

if (!file.exists(z_snp_outfile)){
  # read the data using the VariantAnnotation package
  z_snp <- VariantAnnotation::readVcf(gwas_file)
  z_snp <- as.data.frame(gwasvcf::vcf_to_tibble(z_snp))

  # compute the z-scores
  z_snp$Z <- z_snp$ES/z_snp$SE

  # collect sample size (most frequent sample size for all variants)
  gwas_n <- as.numeric(names(sort(table(z_snp$SS),decreasing=TRUE)[1]))
  cat("gwas_n =", gwas_n, "\n")

  # select example chromosome
  loginfo("Number of SNPs in z_snp: %d", nrow(z_snp))
  z_snp <- z_snp[z_snp$seqnames == example_chrom, ]
  loginfo("Number of SNPs in example chr: %d", nrow(z_snp))

  # subset the columns and format the column names
  z_snp <- z_snp[,c("rsid", "ALT", "REF", "Z")]
  colnames(z_snp) <- c("id", "A1", "A2", "z")

  saveRDS(z_snp, file=z_snp_outfile)
}

## Preprocess z-scores
processed_z_snp_file <- file.path(outputdir, paste0(outname, ".preprocessed.z_snp.RDS"))
if (file.exists(processed_z_snp_file)){
  cat(sprintf("Load preprocessed z_snp: %s \n", processed_z_snp_file))
  z_snp <- readRDS(processed_z_snp_file)
}else{
  z_snp <- readRDS(z_snp_outfile)
  runtime <- system.time({
    z_snp <- preprocess_z_snp(z_snp, snp_info)
  })
  saveRDS(z_snp, processed_z_snp_file)
  loginfo("Preprocessing GWAS z-scores took %0.2f minutes\n",runtime["elapsed"]/60)
}

z_snp_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".preprocessed.z_snp.RDS")))
all.equal(z_snp, z_snp_multigroup)

##### Preprocess weights #####
processed_weight_file <- file.path(outputdir, paste0(outname, ".preprocessed.weights.RDS"))
if (file.exists(processed_weight_file)){
  loginfo("Load preprocessed weight: %s\n", processed_weight_file)
  weights <- readRDS(processed_weight_file)
}else{
  runtime <- system.time({
    weights <- preprocess_weights(weight_file,
                                  region_info = region_info,
                                  gwas_snp_ids = z_snp$id,
                                  snp_info = snp_info,
                                  weight_format = "PredictDB",
                                  ncore = ncore)
  })

  loginfo("Preprocessing weights took %0.2f minutes\n",runtime["elapsed"]/60)
  saveRDS(weights, file = processed_weight_file)
}

weights_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".preprocessed.weights.RDS")))
for (i in 1:length(weights_multigroup)){
  weights_multigroup[[i]]$type <- NULL
  weights_multigroup[[i]]$context <- NULL
}

all.equal(weights, weights_multigroup)

##### Compute gene z-scores #####
cat("##### Computing gene z-scores ##### \n")
gene_z_file <- file.path(outputdir, paste0(outname, ".z_gene.RDS"))
if( file.exists(gene_z_file) ){
  loginfo("Load gene z-scores from %s \n", gene_z_file)
  z_gene <- readRDS(gene_z_file)
}else{
  runtime <- system.time({
    z_gene <- compute_gene_z(z_snp, weights,
                             ncore = ncore,
                             logfile = file.path(outputdir, paste0(outname, ".compute_gene_z.log")))
  })
  saveRDS(z_gene, file = gene_z_file)
  loginfo("Computing gene z-scores took %0.2f minutes\n",runtime["elapsed"]/60)
}

z_gene_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".z_gene.RDS")))
z_gene_multigroup[, c("type", "context", "group")] <- NULL
all.equal(z_gene, z_gene_multigroup)

## Running cTWAS main function
runtime <- system.time({
  ctwas_res <- ctwas_sumstats(z_snp,
                              weights,
                              region_info,
                              snp_info,
                              LD_info,
                              thin = 0.1,
                              maxSNP = 20000,
                              min_nonSNP_PIP = 0.5,
                              ncore = ncore,
                              outputdir = outputdir,
                              outname = paste0(outname, ".ctwas_sumstats"),
                              save_cor = TRUE,
                              cor_dir = file.path(outputdir, "cor_matrix"),
                              logfile = file.path(outputdir, paste0(outname, ".ctwas_sumstats.log")),
                              verbose = FALSE)
})
saveRDS(ctwas_res, file.path(outputdir, paste0(outname, ".ctwas_sumstats_res.RDS")))
loginfo("Running cTWAS main function took %0.2f minutes\n",runtime["elapsed"]/60)

ctwas_res <- readRDS(file.path(outputdir, paste0(outname, ".ctwas_sumstats_res.RDS")))

##### Assemble region_data #####
region_data_file <- file.path(outputdir, paste0(outname, ".region_data.thin", thin, ".RDS"))
boundary_gene_file <- file.path(outputdir, paste0(outname, ".boundary_genes.RDS"))
if (file.exists(region_data_file)) {
  region_data <- readRDS(region_data_file)
  boundary_genes <- readRDS(boundary_gene_file)
} else{
  runtime <- system.time({
    res <- assemble_region_data(region_info,
                                z_snp,
                                z_gene,
                                weights,
                                snp_info,
                                maxSNP = max_snp_region,
                                trim_by = "random",
                                thin = thin,
                                adjust_boundary_genes = TRUE,
                                ncore = ncore)
  })
  loginfo("Assembling region_data took %0.2f minutes\n",runtime["elapsed"]/60)
  region_data <- res$region_data
  boundary_genes <- res$boundary_genes
  saveRDS(region_data, region_data_file)
  saveRDS(boundary_genes, boundary_gene_file)
}

all.equal(ctwas_res$region_data, region_data)
all.equal(ctwas_res$boundary_genes, boundary_genes)

region_data_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".region_data.thin", thin, ".RDS")))
boundary_genes_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".boundary_genes.RDS")))
for (i in 1:length(region_data_multigroup)){
  region_data_multigroup[[i]]$g_type <- NULL
  region_data_multigroup[[i]]$g_context <- NULL
  region_data_multigroup[[i]]$g_group <- NULL
  region_data_multigroup[[i]]$gs_type <- NULL
  region_data_multigroup[[i]]$gs_context <- NULL
  region_data_multigroup[[i]]$gs_group[region_data_multigroup[[i]]$gs_group!="SNP"] <- "gene"
}

all.equal(region_data, region_data_multigroup)
all.equal(boundary_genes, boundary_genes_multigroup)

##### Estimate parameters #####
cat("##### Estimating parameters ##### \n")
# est_param
param_file <- file.path(outputdir, paste0(outname, ".param.RDS"))
if (file.exists(param_file)) {
  param <- readRDS(param_file)
} else{
  runtime <- system.time({
    param <- est_param(region_data,
                       ncore = ncore,
                       verbose = FALSE,
                       logfile = file.path(outputdir, paste0(outname, ".est_param.log")))
  })
  saveRDS(param, param_file)
  loginfo("Parameter estimation took %0.2f minutes\n",runtime["elapsed"]/60)
}
group_prior <- param$group_prior
group_prior_var <- param$group_prior_var
group_size <- param$group_size

all.equal(ctwas_res$param, param)

param_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".param.RDS")))
names(param_multigroup$group_prior) <- c("gene", "SNP")
names(param_multigroup$group_prior_var) <- c("gene", "SNP")
rownames(param_multigroup$group_prior_iters) <- c("gene", "SNP")
rownames(param_multigroup$group_prior_var_iters) <- c("gene", "SNP")
param_multigroup$group_prior_var_structure <- NULL
names(param_multigroup$group_size) <- c("gene", "SNP")
all.equal(param_multigroup, param)

##### Assess parameter estimates #####
ctwas_parameters <- summarize_param(param, gwas_n)
saveRDS(ctwas_parameters, paste0(outputdir, "/", outname, ".ctwas_parameters.RDS"))
ctwas_parameters

ctwas_parameters_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".ctwas_parameters.RDS")))
all.equal(ctwas_parameters, ctwas_parameters_multigroup)

##### Screen regions #####
screen_regions_file <- file.path(outputdir, paste0(outname, ".screened_region_data.RDS"))
if (file.exists(screen_regions_file)) {
  screened_region_data <- readRDS(screen_regions_file)
} else{
  runtime <- system.time({
    screen_regions_res <- screen_regions(region_data,
                                         use_LD = TRUE,
                                         LD_info = LD_info,
                                         snp_info = snp_info,
                                         weights = weights,
                                         group_prior = group_prior,
                                         group_prior_var = group_prior_var,
                                         ncore = ncore,
                                         verbose = FALSE,
                                         logfile = file.path(outputdir, paste0(outname, ".screen_regions.log")))
  })
  loginfo("Screen regions took %0.2f minutes\n",runtime["elapsed"]/60)
  screened_region_data <- screen_regions_res$screened_region_data
  L <- screen_regions_res$L

  # Expand screened region_data with all SNPs in the regions
  if (thin < 1){
    screened_region_data <- expand_region_data(screened_region_data,
                                               snp_info,
                                               z_snp,
                                               z_gene,
                                               trim_by = "z",
                                               maxSNP = max_snp_region,
                                               ncore = ncore)
  }
  saveRDS(screened_region_data, screen_regions_file)
}

all.equal(ctwas_res$screened_region_data, screened_region_data)

screened_region_data_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".screened_region_data.RDS")))
for (i in 1:length(screened_region_data_multigroup)){
  screened_region_data_multigroup[[i]]$g_type <- NULL
  screened_region_data_multigroup[[i]]$g_context <- NULL
  screened_region_data_multigroup[[i]]$g_group <- NULL
  screened_region_data_multigroup[[i]]$gs_type <- NULL
  screened_region_data_multigroup[[i]]$gs_context <- NULL
  screened_region_data_multigroup[[i]]$gs_group[screened_region_data_multigroup[[i]]$gs_group!="SNP"] <- "gene"
}
all.equal(screened_region_data, screened_region_data_multigroup)

##### Finemapping #####
cat("##### Finemapping screened regions ##### \n")

finemap_regions_file <- file.path(outputdir, paste0(outname, ".finemap_regions_res.RDS"))
if (file.exists(finemap_regions_file)) {
  finemap_res <- readRDS(finemap_regions_file)
} else{
  runtime <- system.time({
    finemap_res <- finemap_regions(screened_region_data,
                                   use_LD = TRUE,
                                   LD_info = LD_info,
                                   snp_info = snp_info,
                                   weights = weights,
                                   L = L,
                                   group_prior = group_prior,
                                   group_prior_var = group_prior_var,
                                   ncore = 1,
                                   save_cor = TRUE,
                                   cor_dir = file.path(outputdir, "cor_matrix"),
                                   verbose = TRUE)
  })
  saveRDS(finemap_res, finemap_regions_file)
  loginfo("Finemapping took %0.2f minutes\n",runtime["elapsed"]/60)
}

all.equal(finemap_res, ctwas_res$finemap_res)

finemap_res_multigroup <- readRDS(file.path(multigroup_outputdir, paste0(outname, ".finemap_regions_res.RDS")))
finemap_res_multigroup[,c("type", "context")] <- NULL
finemap_res_multigroup$group[finemap_res_multigroup$group != "SNP"] <- "gene"
all.equal(finemap_res, finemap_res_multigroup)

