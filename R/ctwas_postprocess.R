# get regions with problematic high PIP SNPs or genes
get_problematic_highpip_regions <- function(outputdir, outname, weight, problematic_snps, pip_thresh = 0.5){

  loginfo('Get regions with problematic SNPs')
  loginfo('Number of problematic SNPs: %d', length(problematic_snps))

  if (length(problematic_snps) > 0) {
    # read the PredictDB weights
    stopifnot(file.exists(weight))
    sqlite <- RSQLite::dbDriver("SQLite")
    db <- RSQLite::dbConnect(sqlite, weight)
    query <- function(...) RSQLite::dbGetQuery(db, ...)
    weight_table <- query("select * from weights")
    extra_table <- query("select * from extra")
    # load gene information from PredictDB weights
    gene_info <- query("select gene, genename, gene_type from extra")
    RSQLite::dbDisconnect(db)

    # load cTWAS result
    ctwas_res <- as.data.frame(data.table::fread(file.path(outputdir, paste0(outname, ".susieIrss.txt")), header=T))

    # add gene names to cTWAS results and weight table
    ctwas_res$genename[ctwas_res$type=="gene"] <- gene_info$genename[match(ctwas_res$id[ctwas_res$type=="gene"], gene_info$gene)]

    # find regions with high PIP variants or genes
    ctwas_highpip_res <- ctwas_res[ctwas_res$susie_pip > pip_thresh, ]
    ctwas_highpip_snp_res <- ctwas_highpip_res[ctwas_highpip_res$type == "SNP", ]
    ctwas_highpip_gene_res <- ctwas_highpip_res[ctwas_highpip_res$type == "gene", ]
    ctwas_highpip_gene_weight_table <- weight_table[weight_table$gene %in% ctwas_highpip_gene_res$id, ]

    # find regions with high PIP variants (PIP > 0.5) that are problematic;
    # or high PIP genes with problematic variants in its weights
    problematic_highpip_snps <- intersect(ctwas_highpip_snp_res$id, problematic_snps)
    loginfo('Number of problematic high PIP SNPs: %d', length(problematic_highpip_snps))
    problematic_highpip_genes <- ctwas_highpip_gene_weight_table$gene[which(ctwas_highpip_gene_weight_table$rsid %in% problematic_snps)]
    loginfo('Number of problematic high PIP genes: %d', length(problematic_highpip_genes))

    # get problematic high PIP regions
    ctwas_problematic_res <- ctwas_highpip_res[ctwas_highpip_res$id %in% c(problematic_highpip_snps, problematic_highpip_genes),]
    problematic_regions <- unique(ctwas_problematic_res[,c("region_tag1", "region_tag2")])
    problematic_regions <- problematic_regions[order(problematic_regions$region_tag1,problematic_regions$region_tag2), ]
    problematic_region_tags <- paste0(problematic_regions$region_tag1, "_", problematic_regions$region_tag2)
    loginfo('Number of problematic high PIP regions: %d', length(problematic_region_tags))
  }else{
    loginfo('No problematic SNPs')
    problematic_region_tags <- NULL
  }

  # return a data frame of the problematic region tags
  return(problematic_region_tags)
}

# rerun finemapping with L = 1 for problematic regions
rerun_ctwas_finemap_regions_L1_rss <- function(z_snp,
                                               z_gene,
                                               weight,
                                               problematic_snps = NULL,
                                               rerun_region_tags = NULL,
                                               pip_thresh = 0.5,
                                               ctwas_outputdir,
                                               ctwas_outname,
                                               ld_R_dir = NULL,
                                               outputdir = NULL,
                                               outname = NULL,
                                               group_prior = NULL,
                                               group_prior_var = NULL,
                                               group_prior_var_structure = c("independent","shared_all","shared+snps","inv_gamma","shared_QTLtype"),
                                               inv_gamma_shape=1,
                                               inv_gamma_rate=0,
                                               use_null_weight = T,
                                               coverage = 0.95,
                                               min_abs_corr = 0.5,
                                               ncore = 1,
                                               logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file= logfile, level='DEBUG')
  }

  # get regions with problematic high PIP SNPs or genes
  if (is.null(rerun_region_tags) && length(problematic_snps) > 0) {
    rerun_region_tags <- get_problematic_highpip_regions(outputdir = ctwas_outputdir,
                                                         outname = ctwas_outname,
                                                         weight = weight,
                                                         problematic_snps = problematic_snps,
                                                         pip_thresh = pip_thresh)
  }

  if (length(rerun_region_tags) > 0) {

    if (is.null(outputdir)) {
      outputdir <- paste0(ctwas_outputdir, "/rerun_regions/")
    }
    if (is.null(outname)) {
      outname <- ctwas_outname
    }
    dir.create(outputdir, showWarnings=F, recursive = T)

    # Load cTWAS result
    ctwas_res <- as.data.frame(data.table::fread(file.path(ctwas_outputdir, paste0(ctwas_outname, ".susieIrss.txt")), header=T))

    # Load estimated parameters
    if (is.null(group_prior) || is.null(group_prior_var)) {
      loginfo("Load estimated parameters")
      load(file.path(ctwas_outputdir, paste0(ctwas_outname, ".s2.susieIrssres.Rd")))
      group_prior <- group_prior_rec[,ncol(group_prior_rec)]
      group_prior_var <- group_prior_var_rec[,ncol(group_prior_var_rec)]
      rm(group_prior_rec, group_prior_var_rec)
    }
    group_prior_var_structure <- match.arg(group_prior_var_structure)

    # Load regionlist
    regionlist <- readRDS(paste0(ctwas_outputdir, "/", ctwas_outname, ".regionlist.RDS"))
    regs <- data.table::fread(paste0(ctwas_outputdir,"/", ctwas_outname, ".regions.txt"))

    # select and assemble regionlist for rerunning finemapping
    subset_regionlist_res <- subset_regionlist(regionlist, rerun_region_tags)
    rerun_regionlist <- subset_regionlist_res$regionlist
    rerun_regs <- subset_regionlist_res$regs

    saveRDS(rerun_regionlist, paste0(outputdir, "/",outname, ".rerun.regionlist.RDS"))
    write.table(rerun_regs, paste0(outputdir,"/", outname, ".rerun.regions.txt"),
                row.names=F, col.names=T, sep="\t", quote = F)

    # position information for imputed genes
    ld_exprvarfs <- paste0(ctwas_outputdir, "/", ctwas_outname, "_chr", 1:22, ".exprvar")

    if(is.null(ld_R_dir)) {
      ld_R_dir <- dirname(regionlist[[1]][[1]]$regRDS)[1]
    }
    loginfo("ld_R_dir: %s", ld_R_dir)

    loginfo('Rerun finemmapping for %d regions', nrow(rerun_regs))
    loginfo('Output directory: %s', outputdir)

    # Rerun finemapping with L = 1
    ctwas_finemap_rss(z_gene = z_gene,
                      z_snp = z_snp,
                      ld_exprvarfs = ld_exprvarfs,
                      ld_R_dir = ld_R_dir,
                      reuse_regionlist = T,
                      regionlist_custom = rerun_regionlist,
                      outputdir = outputdir,
                      outname = paste0(outname, ".L1"),
                      L = 1,
                      group_prior = group_prior,
                      group_prior_var = group_prior_var,
                      group_prior_var_structure = group_prior_var_structure,
                      use_null_weight = use_null_weight,
                      coverage = coverage,
                      min_abs_corr = min_abs_corr,
                      inv_gamma_shape=inv_gamma_shape,
                      inv_gamma_rate=inv_gamma_rate,
                      ncore = ncore)

  }

}

# select and assemble a subset of regionlist by region_tags
subset_regionlist <- function(regionlist, region_tags){

  # subset regionlist
  region_tags_df <- data.frame(
    b = sapply(strsplit(region_tags, "_"), "[[", 1),
    rn = sapply(strsplit(region_tags, "_"), "[[", 2))

  regionlist_subset <- vector("list", length = 22)
  for (b in 1:22) {
    rn <- region_tags_df[region_tags_df$b == b, "rn"]
    regionlist_subset[[b]] <- regionlist[[b]][as.character(rn)]
  }
  region_subset_chrs <- sort(unique(as.integer(region_tags_df$b)))

  if(length(regionlist_subset) > 0){
    # coordinates of the selected regions
    temp_regs <- lapply(1:22, function(x) cbind(x,
                                                unlist(lapply(regionlist_subset[[x]], "[[", "start")),
                                                unlist(lapply(regionlist_subset[[x]], "[[", "stop"))))

    regs <- do.call(rbind, lapply(temp_regs, function(x) if (ncol(x) == 3){x}))
    loginfo("select %d regions from the regionlist", nrow(regs))
  }else{
    loginfo("No regions selected.")
    regionlist_subset <- NULL
    regs <- NULL
  }

  return(list(regionlist = regionlist_subset,
              regs = regs))
}
