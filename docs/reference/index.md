# Package index

## All functions

<!-- end list -->

  - `add_pos()` : Add SNP and gene positions to cTWAS finemapping
    result.
  - `anno_finemap_res()` : Map finemapping result of molecular traits to
    genes.
  - `anno_susie_alpha_res()` : Map molecular traits to genes in susie
    alpha result.
  - `assemble_region_data()` : Assembles data for all the regions
  - `check_n_snps()` : Check the numbers of SNPs in snp\_map, z\_snp and
    weights
  - `combine_gene_pips()` : Combines gene PIPs by context, type or
    group.
  - `compute_enrichment_test()` : Computes enrichment (log-scale),
    standard error and p-value
  - `compute_gene_z()` : Computes z-scores of molecular traits
  - `compute_region_cor()` : Computes correlation matrices for a single
    region
  - `compute_region_nonSNP_PIPs()` : Computes non-SNP PIPs for all
    regions from finemapping result
  - `compute_weight_LD_from_ref()` : Computes LD for weight variants
    using reference LD
  - `convert_fusion_to_predictdb()` : Converts fusion weights to
    predictDB format
  - `convert_geno_to_LD_matrix()` : Converts PLINK genotype data to LD
    matrices and SNP info files, saves LD matrices as .RDS files and SNP
    info as .Rvar files
  - `convert_predictdb_weights_varIDs()` : Convert rsIDs and varIDs in
    predictDB weights to the reference variant format.
  - `convert_to_ukb_varIDs()` : Convert variant IDs from Open GWAS
    format or PredictDB weight format
  - `create_predictdb_from_QTLs()` : Creates weight files in PredictDB
    format from QTL data
  - `create_snp_LD_map()` : Map SNPs to regions using region meta table.
  - `create_snp_map()` : Map SNPs to regions using all the variants in
    the LD reference.
  - `ctwas_sumstats()` : cTWAS analysis using summary statistics
  - `ctwas_sumstats_noLD()` : cTWAS analysis using summary statistics
    with "no LD" version
  - `diagnose_LD_mismatch_susie()` : Diagnose LD mismatch using SuSiE
    RSS
  - `est_param()` : Estimates cTWAS parameters using EM with cTWAS SER
    model
  - `estimate_region_L()` : Estimate L for all regions by running
    finemapping with uniform prior
  - `expand_region_data()` : Expands region\_data with all SNPs
  - `filter_z_gene_by_group_size()` : Filter z\_gene by group size
  - `finemap_regions()` : Runs cTWAS fine-mapping for regions
  - `finemap_regions_noLD()` : Runs cTWAS fine-mapping for regions
    without LD (L = 1)
  - `get_boundary_genes()` : Get cross-boundary genes.
  - `get_gene_annot_from_ens_db()` : Get gene annotation table from
    Ensembl database
  - `get_gene_info()` : Get gene info from weights.
  - `get_molecular_ids()` : Get original molecular IDs
  - `get_predictdb_genome_build()` : Get genome build of PredictDB
    weight
  - `get_problematic_genes()` : Gets problematic genes from problematic
    SNPs
  - `get_region_cor()` : Gets correlation matrices for a single region.
  - `load_LD()` : Load LD matrix
  - `load_fusion_weights()` : Loads weights in FUSION format
  - `load_predictdb_weights()` : Loads weights in PredictDB format
  - `load_region_cor()` : Loads precomputed correlation matrices for a
    single region.
  - `load_weights()` : Loads weights in PredictDB or FUSION format
  - `make_convergence_plots()` : Make convergence plots for the
    estimated parameters
  - `make_locusplot()` : Plots the cTWAS result for a single locus
  - `map_gene_regions()` : Map regions for each gene.
  - `merge_region_data()` : Merges region data for cross-boundary genes
  - `merge_region_data_noLD()` : Merges region data for cross-boundary
    genes without using LD
  - `postprocess_LD_mismatch()` : Runs cTWAS post-processing procedure
    for merging regions
  - `postprocess_region_merging()` : Runs cTWAS post-processing
    procedure for region merging
  - `postprocess_region_merging_noLD()` : Runs cTWAS post-processing
    procedure for region merging without LD
  - `preprocess_weights()` : Preprocess PredictDB/FUSION weights and
    harmonize with LD reference
  - `preprocess_z_snp()` : Preprocess GWAS z-scores, harmonize GWAS
    z-scores with LD reference
  - `read_gwas()` : Read GWAS summary statistics
  - `read_snp_info_file()` : Read a single SNP info file as a data frame
  - `read_snp_info_files()` : Read multiple single SNP info files as a
    data frame
  - `screen_regions()` : Screens regions with strong signals
  - `select_boundary_genes()` : Gets boundary genes and selects high PIP
    boundary genes
  - `subset_weights()` : Get a subset of weights, by group, context or
    type.
  - `summarize_param()` : Summarizes estimated parameters
  - `summarize_region_signals()` : Gets a basic summary of region
    signals
  - `trim_weights()` : Trim weights and filter SNPs by LD reference and
    GWAS SNPs.
  - `update_finemap_res()` : Updates cTWAS finemapping result for
    selected regions
  - `update_merged_region_data()` : Updates cTWAS input data with merged
    region data
  - `update_merged_region_data_noLD()` : Updates cTWAS input data
    without LD with merged region data
  - `update_merged_region_finemap_res()` : Updates cTWAS finemapping
    result for merged regions
  - `update_region_z()` : Adds or updates z-scores in region\_data based
    on z\_snp and z\_gene. this will also update sid and gid based on
    z\_snp and z\_gene.
  - `z2p()` : convert z-scores to p-values (two-sided test)
