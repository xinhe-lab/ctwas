#' Causal inference for TWAS using summary statistics
#'
#' @param z_snp A data frame with four columns: "id", "A1", "A2", "z".
#' giving the z scores for snps. "A1" is effect allele. "A2" is the other allele.
#'
#' @param weight a string, pointing to a directory with the fusion/twas format of weights, or a .db file in predictdb format.
#' A vector of multiple sets of weights in PredictDB format can also be specified; genes will have their filename appended
#' to their gene name to ensure IDs are unique.

#' @param region_info a data frame of region definition and associated file names
#'
#' @param outputdir a string, the directory to store output
#'
#' @param outname a string, the output name
#'
#' @param logfile the log file, if NULL will print log info on screen
#'
#' @importFrom logging addHandler loginfo writeToFile
#'
#' @return a list of imputed gene z-scores, estimated parameters,
#' finemapping PIPs and credible sets, and updated region_info
#'
#' @export
#'
ctwas_sumstats <- function(
    z_snp,
    weight,
    region_info,
    outputdir = getwd(),
    outname = "ctwas_sumstats",
    logfile = NULL){

  if (!is.null(logfile)){
    addHandler(writeToFile, file=logfile, level='DEBUG')
  }

  # compute gene z-scores
  z_gene <- compute_gene_z(z_snp, weight, region_info)

  # estimate parameters (including compute_cor)
  param <- est_param(z_snp, z_gene, region_info, outputdir, outname)

  # screen regions
  regionlist <- screen_regions(z_snp, z_gene, param, region_info)

  # fine-map selected regions (including compute_cor)
  fine_map_regions(z_snp, z_gene, param, region_info, regionlist)

  return(list("z_gene" = z_gene,
              "param" = par,
              "ctwas_susie_res" = ctwas_susie_res,
              "region_info" = updated_region_info))

}

