
# Example tables:
# 'extra' table
# gene|genename|gene_type|n.snps.in.model|pred.perf.R2|pred.perf.pval|pred.perf.qval
# ENSG00000000457.13|SCYL3|protein_coding|2|||
# ENSG00000000460.16|C1orf112|protein_coding|2|||
# ENSG00000000938.12|FGR|protein_coding|1|||
# ENSG00000000971.15|CFH|protein_coding|1|||
# ENSG00000001036.13|FUCA2|protein_coding|2|||

# 'weights' table
# gene|rsid|varID|ref_allele|eff_allele|weight
# ENSG00000169583.12|rs141364387|chr9_136998041_C_T_b38|C|T|-0.15870216858425
# ENSG00000107331.16|rs57213430|chr9_137029055_C_CA_b38|C|CA|0.197539699546426
# ENSG00000180549.7|rs10732706|chr9_137032508_C_T_b38|C|T|-0.177871912734104
# ENSG00000107281.9|rs2891950|chr9_137045741_C_G_b38|C|G|0.539136477249632
# ENSG00000107281.9|chr9_137046201_C_A_b38|chr9_137046201_C_A_b38|C|A|0.320361972555521 
main <- function(
    template_db_path = "../input/template_predixcan.db",
    rdata_dir = "../input/wingo_nc_brain_expression_888_weights_raw",
    out_db_path = "../output/wingo_nc_brain_expression_888_predixcan.db",
    expect_hg38_pos_in_snps=F,
    want_model_type = "lasso" # Choose from "blup"  "lasso" "top1"  "enet"  "bslmm"
){
  require(tidyverse)
  require(purrr)
  require(glue)
  
  # Copy the template db to the new db name
  out_db_path <- glue(out_db_path)
  file.copy(template_db_path, out_db_path, overwrite = TRUE)
  
  # Get list of all files in rdata_dir
  rdata_files <- list.files(rdata_dir, full.names = TRUE)
  print(length(rdata_files))
  process_rdata_file <- function(rdata_fpath){
    # Filenames have form: {ensg_no_decimal}.wgt.RDat
    ensg_no_dec <- basename(rdata_fpath) %>% str_replace_all(".wgt.RDat", "")
    load(rdata_fpath)
    
    # Get the best performing model (lowest p-val/highest rsq)
    cv.rsq <- cv.performance[1,]
    want_model_index <- which(names(cv.rsq) == want_model_type)
    want_model_rsq <- cv.rsq[want_model_index]
    
    want_model_weights_df <- as_tibble(wgt.matrix[,want_model_index]) %>%
      mutate(rsid = rownames(wgt.matrix)) %>%
      filter(value != 0) %>%
      left_join(as_tibble(snps) %>% 
                  rename(chr=V1,
                         rsid=V2,
                         cm=V3,
                         pos_hg37=V4,
                         ref_allele=V6, # Taking the first column in bim file analog to be effect allele
                         eff_allele=V5) %>%
                  select(-cm)
      )
    if (nrow(want_model_weights_df) == 0){
      outcome <- "found_all_0_weights_model"
    } else {
      outcome <- glue("found_{want_model_type}_model")
    }
    
    # If the position in the SNPs file is actually hg38, then we can use a "proper" varID
    # so that we can calculate gtex ld if needed
    if (expect_hg38_pos_in_snps==F){
      want_model_weights_insert_df <- want_model_weights_df %>%
        rename(weight = value) %>%
        mutate(varID = glue("{chr}_{pos_hg37}"),
               gene = ensg_no_dec) %>%
        select(all_of(c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")))
    } else { # Else, pos_hg37 is actually pos_hg38 so we can use full SPrediXcan style varID column
      want_model_weights_insert_df <- want_model_weights_df %>%
        rename(weight = value) %>%
        mutate(varID = glue("chr{chr}_{pos_hg37}_{ref_allele}_{eff_allele}_b38"),
               gene = ensg_no_dec) %>%
        select(all_of(c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")))
    }
    
    # Insert into 'extra' table
    # 'extra' table
    # gene|genename|gene_type|n.snps.in.model|pred.perf.R2|pred.perf.pval|pred.perf.qval
    out_db_con <- DBI::dbConnect(RSQLite::SQLite(), dbname = out_db_path)
    append_extra_res <- DBI::dbAppendTable(out_db_con, "extra",
                                           tibble(gene = ensg_no_dec,
                                                  genename = NA,
                                                  gene_type = NA,
                                                  n.snps.in.model = nrow(want_model_weights_df),
                                                  pred.perf.R2 = want_model_rsq,
                                                  pred.perf.pval = NA,
                                                  pred.perf.qval = NA
                                           )
    )
    
    # Insert into 'weights' table
    # 'weights' table
    # gene|rsid|varID|ref_allele|eff_allele|weight
    # ENSG00000000457.13|SCYL3|protein_coding|2|||
    # ENSG00000169583.12|rs141364387|chr9_136998041_C_T_b38|C|T|-0.15870216858425
    append_weights_res <- DBI::dbAppendTable(out_db_con, "weights", want_model_weights_insert_df)
    DBI::dbDisconnect(out_db_con)
    
    # Check that the inserts worked
    if (append_extra_res == 1 & append_weights_res == nrow(want_model_weights_insert_df)){
      print(glue("Inserts worked for {ensg_no_dec}."))
    } else {
      return(tibble(group=ensg_no_dec, outcome="insertion_failure"))
      message(glue("One or both inserts failed for {ensg_no_dec}."))
    }
    return(tibble(group=ensg_no_dec, outcome=outcome))
  }
  
  outcome_df <- rdata_files %>% map_dfr(process_rdata_file)
  
  # Write outcome dataframe
  outcome_df_name <- out_db_path %>%
    str_replace_all("\\.db", ".tsv")
  
  write_delim(outcome_df, outcome_df_name, delim="\t")
}

# main <- function(
    #   template_db_path = "../input/template_predixcan.db",
#   rdata_dir = "../input/wingo_nc_brain_expression_888_weights_raw",
#   out_db_path = "../output/wingo_nc_brain_expression_888_predixcan.db"
#   ){
if (!interactive()) {
  print("Sourced non-interactively.")
  # main(rdata_dir = "../input/nilanchatterjee_pwas_Plasma_Protein_weights_EA",
  #      out_db_path = "../output/nilanchatterjee_pwas_EA_predixcan.db",
  #      expect_hg38_pos_in_snps=T
  #      )
  
  # main(rdata_dir = "../input/Wingo_NC_2022_pqtl_fusion_weights",
  #      out_db_path = "../output/Wingo_NC_2022_pqtl_fusion_{want_model_type}_models_predixcan.db")
} else {
  print("Sourced interactively.")
  main(rdata_dir = "../input/Wingo_NC_2022_pqtl_fusion_weights",
       out_db_path = "../output/Wingo_NC_2022_pqtl_fusion_{want_model_type}_models_predixcan.db")
}
