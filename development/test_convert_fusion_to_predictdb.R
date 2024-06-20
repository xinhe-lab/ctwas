
source("/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/fusiontodb_julian_JG.R")
fusion_to_predictDB(
  rdata_dir = "/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/fusion_weights",
  template_db_path = "/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/template_predixcan.db",
  out_db_path = "/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/test.db",
  expect_hg38_pos_in_snps=FALSE,
  want_model_type = "lasso")

res <- convert_fusion_to_predictdb(
  weight_path = "/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/fusion_weights",
  fusion_method = "lasso",
  fusion_genome_version = "b37",
  outputdir = "/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/",
  outname = "test2")

res <- load_predictdb_weights("/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/test.db")

res2 <- load_predictdb_weights("/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/test2.db")

res <- convert_fusion_to_predictdb(
  weight_file = "/project/mstephens/causalTWAS/apa_models/Heart_Atrial_Appendage/Heart_Atrial_Appendage",
  fusion_method = "bestR2",
  fusion_genome_version = "b38",
  outputdir = "/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/",
  outname = "test_Heart_Atrial_Appendage")

res <- load_predictdb_weights("/project2/xinhe/shared_data/ctwas_tutorial/scripts/fusion_to_db/test_Heart_Atrial_Appendage.db")
