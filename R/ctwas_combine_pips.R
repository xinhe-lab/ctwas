
#' @title Combines gene PIPs by context, type or group.
#'
#' @param susie_alpha_res a data frame of annotated susie alpha result.
#'
#' @param group_by column name to group genes by.
#'
#' @param by option to combine PIPs by: "context" (default), "type", or "group".
#'
#' @param method method to combine PIPs of molecular traits targeting the same gene.
#' options:
#' "combine_cs" (default):
#' first sums PIPs of molecular traits of a genes in each credible set,
#' and then combine PIPs using the following formula:
#' \eqn{1 - \prod_k (1 - \text{PIP}_k)},
#' where \eqn{\text{PIP}_k} is the summed PIP of the \eqn{k}-th credible set of a gene.
#' This is the default option for combining PIPs from fine-mapping with LD.
#' "sum": sum over PIPs of all molecular traits for the same gene.
#' This summation is the expected number of causal molecular traits in this gene,
#' and could be higher than 1.
#' We will use this option for combining PIPs from fine-mapping without LD.
#'
#' @param filter_cs If TRUE, limits gene results to credible sets (CS).
#'
#' @param keep_alpha_in_cs_only If TRUE, only keep single effects (alpha) in credible sets
#' when calculating combined PIP. This is similar to \code{prune_by_cs} in susie.
#'
#' @param include_cs_id If TRUE, include credible set IDs of the genes in the output
#'
#' @param include_set_id If TRUE, include susie set IDs of the genes in the output
#'
#' @param missing_value set missing value as (default: NA)
#'
#' @return a data frame of combined gene PIPs for each context, type or group
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#'
#' @export
combine_gene_pips <- function(susie_alpha_res,
                              group_by = "molecular_id",
                              by = c("context", "type", "group"),
                              method = c("combine_cs", "sum"),
                              filter_cs = TRUE,
                              keep_alpha_in_cs_only = FALSE,
                              include_cs_id = TRUE,
                              include_set_id = FALSE,
                              missing_value = NA){

  by <- match.arg(by)
  method <- match.arg(method)

  # Check to see if gene_name and gene_type are already in susie_alpha_res
  if (!(group_by %in% colnames(susie_alpha_res))){
    stop(paste("Cannot find the column", group_by, "in susie_alpha_res!"))
  }

  # use sum option if there is no CS information
  if (is.null(susie_alpha_res$cs)){
    filter_cs <- FALSE
    include_cs_id <- FALSE
    method <- "sum"
  }

  # Only keep gene results
  susie_alpha_res <- susie_alpha_res[susie_alpha_res$group!="SNP",,drop=FALSE]

  # combine PIPs across all contexts or types
  combined_gene_pips <- compute_combined_pips(susie_alpha_res,
                                              group_by = group_by,
                                              method = method,
                                              filter_cs = filter_cs,
                                              keep_alpha_in_cs_only = keep_alpha_in_cs_only,
                                              include_set_id = include_set_id,
                                              include_cs_id = include_cs_id)
  if (include_set_id) {
    colnames(combined_gene_pips)[colnames(combined_gene_pips) == "set_id"] <- "combined_set_id"
  }

  if (include_cs_id) {
    colnames(combined_gene_pips)[colnames(combined_gene_pips) == "cs_id"] <- "combined_cs_id"
  }

  if (by == "context"){
    # combine PIPs for each context
    contexts <- unique(susie_alpha_res$context)
    for (context in contexts){
      tmp_susie_alpha_res <- susie_alpha_res[susie_alpha_res$context==context,,drop=FALSE]
      tmp_combined_gene_pips <- compute_combined_pips(tmp_susie_alpha_res,
                                                      group_by = group_by,
                                                      method = method,
                                                      filter_cs = filter_cs,
                                                      keep_alpha_in_cs_only = keep_alpha_in_cs_only,
                                                      include_set_id = include_set_id,
                                                      include_cs_id = include_cs_id)

      colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "combined_pip"] <-
        paste0(context, "_pip")

      if (include_set_id) {
        colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "set_id"] <-
          paste0(context, "_set_id")
      }

      if (include_cs_id) {
        colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "cs_id"] <-
          paste0(context, "_cs_id")
      }

      combined_gene_pips <- combined_gene_pips %>%
        left_join(tmp_combined_gene_pips, by = group_by)
    }
  } else if (by == "type") {
    # combine PIPs for each type
    types <- unique(susie_alpha_res$type)
    for (type in types){
      tmp_susie_alpha_res <- susie_alpha_res[susie_alpha_res$type==type,,drop=FALSE]
      tmp_combined_gene_pips <- compute_combined_pips(tmp_susie_alpha_res,
                                                      group_by = group_by,
                                                      method = method,
                                                      filter_cs = filter_cs,
                                                      keep_alpha_in_cs_only = keep_alpha_in_cs_only,
                                                      include_set_id = include_set_id,
                                                      include_cs_id = include_cs_id)

      colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "combined_pip"] <-
        paste0(type, "_pip")

      if (include_set_id) {
        colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "set_id"] <-
          paste0(type, "_set_id")
      }

      if (include_cs_id) {
        colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "cs_id"] <-
          paste0(type, "_cs_id")
      }

      combined_gene_pips <- combined_gene_pips %>%
        left_join(tmp_combined_gene_pips, by = group_by)
    }
  } else if (by == "group") {
    # combine PIPs for each group
    groups <- unique(susie_alpha_res$group)
    for (group in groups){
      tmp_susie_alpha_res <- susie_alpha_res[susie_alpha_res$group==group,,drop=FALSE]
      tmp_combined_gene_pips <- compute_combined_pips(tmp_susie_alpha_res,
                                                      group_by = group_by,
                                                      method = method,
                                                      filter_cs = filter_cs,
                                                      keep_alpha_in_cs_only = keep_alpha_in_cs_only,
                                                      include_set_id = include_set_id,
                                                      include_cs_id = include_cs_id)

      colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "combined_pip"] <-
        paste0(group, "_pip")

      if (include_set_id) {
        colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "set_id"] <-
          paste0(group, "_set_id")
      }

      if (include_cs_id) {
        colnames(tmp_combined_gene_pips)[colnames(tmp_combined_gene_pips) == "cs_id"] <-
          paste0(group, "_cs_id")
      }

      combined_gene_pips <- combined_gene_pips %>%
        left_join(tmp_combined_gene_pips, by = group_by)
    }
  }

  if (!is.na(missing_value)) {
    combined_gene_pips[is.na(combined_gene_pips)] <- missing_value
  }

  # order by combined PIP
  combined_gene_pips <- combined_gene_pips[order(-combined_gene_pips$combined_pip),]
  rownames(combined_gene_pips) <- NULL

  # new_colnames <- c(setdiff(colnames(combined_gene_pips), "combined_pip"), "combined_pip")
  # combined_gene_pips <- combined_gene_pips[, new_colnames]

  return(combined_gene_pips)
}


# Computes combined gene PIPs.
# "combine_cs" (default): first sum alpha of molecular traits for the same gene in the same CS,
# then apply the multiplication formula across CS.
# "sum" sum over alpha of all molecular traits for the same gene;
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join group_by summarise
#'
compute_combined_pips <- function(susie_alpha_res,
                                  group_by = "molecular_id",
                                  method = c("combine_cs", "sum"),
                                  filter_cs = TRUE,
                                  keep_alpha_in_cs_only = FALSE,
                                  include_cs_id = TRUE,
                                  include_set_id = FALSE){

  method <- match.arg(method)

  # Only keep gene results
  susie_alpha_res <- susie_alpha_res[susie_alpha_res$group!="SNP",,drop=FALSE]

  susie_alpha_res$set_id <- paste0(susie_alpha_res$region_id, ".", susie_alpha_res$susie_set)
  susie_alpha_res$cs_id <- paste0(susie_alpha_res$region_id, ".", susie_alpha_res$cs)
  susie_alpha_res$cs_id[which(is.na(susie_alpha_res$cs))] <- ""

  if (filter_cs) {
    # limit to genes in credible sets
    susie_alpha_res <- susie_alpha_res[!is.na(susie_alpha_res$cs),,drop=FALSE]
  }

  if (keep_alpha_in_cs_only) {
    # only keep alpha in credible sets
    susie_alpha_res <- susie_alpha_res[susie_alpha_res$in_cs == TRUE,,drop=FALSE]
  }

  if (method == "sum") {
    combined_gene_pips <- susie_alpha_res %>%
      group_by(.data[[group_by]]) %>%
      summarise(set_id = paste(unique(.data$set_id), collapse = ","),
                cs_id = paste(unique(.data$cs_id[.data$cs_id!=""]), collapse = ","),
                combined_pip = sum(.data$susie_alpha))

  } else if (method == "combine_cs") {
    # for each gene, first sum PIPs within the same sets
    summed_gene_pips <- susie_alpha_res %>%
      group_by(.data[[group_by]], .data$set_id) %>%
      summarise(summed_pip = sum(.data$susie_alpha),
                cs_id = paste(unique(.data$cs_id[.data$cs_id!=""]), collapse = ","))

    # then combine PIPs using the multiplication formula across sets
    combined_gene_pips <- summed_gene_pips %>%
      group_by(.data[[group_by]]) %>%
      summarise(set_id = paste(unique(.data$set_id), collapse = ","),
                cs_id = paste(unique(.data$cs_id[.data$cs_id!=""]), collapse = ","),
                combined_pip = combine_pips_fun(.data$summed_pip))
  }

  combined_gene_pips <- as.data.frame(combined_gene_pips)

  if (!include_set_id) {
    combined_gene_pips <- combined_gene_pips[, colnames(combined_gene_pips)!="set_id"]
  }

  if (!include_cs_id) {
    combined_gene_pips <- combined_gene_pips[, colnames(combined_gene_pips)!="cs_id"]
  }

  return(combined_gene_pips)
}


# Compute combined gene PIP using the multiplication formula
# combined gene PIP = 1 - \prod_k (1 - PIP_k).
# PIP_k is the PIP of the k-th molecular trait of a gene.
combine_pips_fun <- function(pips){
  return(1 - prod(1 - pips))
}

