
#' Plot the cTWAS result for a single locus
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param region_id region ID to be plotted
#'
#' @param weights a list of weights
#'
#' @param ens_db Ensembl database
#'
#' @param use_LD TRUE/FALSE. If TRUE, use LD for finemapping. Otherwise, use "no-LD" version.
#'
#' @param LD_info a list of paths to LD matrices for each of the regions. Required when \code{use_LD = TRUE}.
#'
#' @param snp_info a list of SNP info data frames for LD reference. Required when \code{use_LD = TRUE}.
#'
#' @param save_cor TRUE/FALSE. If TRUE, save correlation (R) matrices to \code{cor_dir}
#'
#' @param cor_dir path of correlation matrices
#'
#' @param region_data a list of region_data used for finemapping
#'
#' @param locus_range a vector of start and end positions to define the locus region to be plotted
#'
#' If NULL, the entire region will be plotted (with 100 bp flanks on both sides).
#'
#' @param focus_gene the gene name to focus the plot. If NULL, choose the gene with the highest PIP
#'
#' @param gene_biotype specified biotypes to be displayed in gene tracks or NULL to display all possible ones.
#' By default, only show protein coding genes.
#'
#' @param point.sizes size values for SNP and non-SNP data points in the scatter plots
#'
#' @param point.alpha alpha values for SNP and non-SNP data points in the scatter plots
#'
#' @param point.shapes shapes values for data points of different types in the scatter plots
#'
#' @param genelabel.size Font size for gene text
#'
#' @param legend.text.size Font size for legend text
#'
#' @param legend.position position to put legends. If "none", no legends will be shown.
#'
#' @param panel.heights Relative heights of the panels.
#'
#' @param highlight_pos If not NULL, highlighst the position in all panels for the genomic position.
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import ggrepel
#' @importFrom locuszoomr locus gg_genetracks
#' @importFrom logging loginfo
#'
#' @export
#'
make_locusplot <- function(finemap_res,
                           region_id,
                           weights,
                           ens_db,
                           use_LD = TRUE,
                           LD_info = NULL,
                           snp_info = NULL,
                           save_cor = FALSE,
                           cor_dir = NULL,
                           region_data = NULL,
                           locus_range = NULL,
                           focus_gene = NULL,
                           gene_biotype = "protein_coding",
                           p_thresh = NULL,
                           pip_thresh = NULL,
                           point.sizes = c(1, 3.5),
                           point.alpha = c(0.4, 0.6),
                           point.shapes = c(21:25),
                           genelabel.size = 2.5,
                           legend.text.size = 10,
                           legend.position = "none",
                           panel.heights = c(4, 4, 0.3, 4),
                           highlight_pos = NULL) {

  if (!is.list(weights)){
    stop("'weights' should be a list.")
  }

  # select finemapping result for the target region
  finemap_region_res <- finemap_res[which(finemap_res$region_id==region_id), ]
  finemap_region_res$p <- (-log(2) - pnorm(abs(finemap_region_res$z), lower.tail=F, log.p=T))/log(10) # add pvalue
  finemap_region_res$object_type <- finemap_region_res$type # adding object_type for size and alpha. SNP and non-SNP
  finemap_region_res$object_type[finemap_region_res$object_type!="SNP"] <- "non-SNP"

  if (!is.null(focus_gene)) {
    focus_gidx <- which(finemap_region_res$gene_name == focus_gene)
    if (length(focus_gidx) == 0){
      stop("can't find focus_gene in the fine-mapping result of the region")
    } else if (length(focus_gidx) > 1){
      finemap_region_gene_res <- finemap_region_res[focus_gidx,]
      focus_gid <- finemap_region_gene_res$id[which.max(finemap_region_gene_res$susie_pip)]
    } else{
      focus_gid <- finemap_region_res$id[focus_gidx]
    }
  } else{
    # if focus_gid not specified, select the top gene
    finemap_region_gene_res <- finemap_region_res[finemap_region_res$type != "SNP",]
    focus_gidx <- which.max(finemap_region_gene_res$susie_pip)
    focus_gene <- finemap_region_gene_res$gene_name[focus_gidx]
    focus_gid <- finemap_region_gene_res$id[focus_gidx]
  }
  loginfo("focus gene: %s", focus_gene)
  loginfo("focus id: %s", focus_gid)

  if (!is.null(locus_range)){
    if (is.null(finemap_region_res$pos)){
      stop("please annotate finemapping result first!")
    }
    idx_in_range <- which(finemap_region_res$pos>=locus_range[1] & finemap_region_res$pos<=locus_range[2])
    finemap_region_res <- finemap_region_res[idx_in_range,, drop=F]
  }else{
    # if locus_range not specified, plot the whole region
    locus_range <- c(min(finemap_region_res$pos)-100, max(finemap_region_res$pos) + 100)
  }
  loginfo("locus_range: %s", locus_range)

  if (use_LD) {
    # add r2 with the focus gene
    # load precomputed correlation matrices if available
    if (!is.null(cor_dir)) {
      if (!dir.exists(cor_dir))
        dir.create(cor_dir, recursive = TRUE)
      R_sg_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp_gene.RDS"))
      R_g_file <- file.path(cor_dir, paste0("region.", region_id, ".R_gene.RDS"))
      R_s_file <- file.path(cor_dir, paste0("region.", region_id, ".R_snp.RDS"))
    }

    cor_files_exist <- isTRUE(!is.null(cor_dir) && file.exists(R_sg_file) && file.exists(R_g_file) && file.exists(R_s_file))

    if (cor_files_exist) {
      # load precomputed correlation matrices
      loginfo("Load correlation matrices for region %s", region_id)
      R_snp_gene <- load_LD(R_sg_file)
      R_gene <- load_LD(R_g_file)
      R_snp <- load_LD(R_s_file)

      if (!is.null(region_data)) {
        gids <- region_data[[region_id]]$gid
        sids <- region_data[[region_id]]$sid
        if (region_data[[region_id]]$thin != 1){
          stop("thin != 1 in the region_data, please expand the region with full SNPs")
        }
        rownames(R_snp) <- sids
        colnames(R_snp) <- sids
        rownames(R_snp_gene) <- sids
        colnames(R_snp_gene) <- gids
        rownames(R_gene) <- gids
        colnames(R_gene) <- gids
      }
    } else {
      # compute correlations
      loginfo("Compute correlation matrices for region %s", region_id)

      if (is.null(region_data)) {
        stop("region_data is needed for computing correlation matrices")
      }
      if (region_data[[region_id]]$thin != 1){
        stop("thin != 1 in the region_data, please expand the region with full SNPs")
      }
      # load LD matrix of the region
      if (is.null(LD_info) || is.null(snp_info)) {
        stop("LD_info and snp_info are required for computing correlation matrices")
      }
      LD_matrix_files <- LD_info[[region_id]]$LD_matrix
      stopifnot(all(file.exists(LD_matrix_files)))
      if (length(LD_matrix_files)==1) {
        R_snp <- load_LD(LD_matrix_files)
      } else {
        R_snp <- lapply(LD_matrix_files, load_LD)
        R_snp <- suppressWarnings(as.matrix(Matrix::bdiag(R_snp)))
      }
      # load SNP info of the region
      snpinfo <- snp_info[[region_id]]

      # compute correlation matrices
      sids <- region_data[[region_id]]$sid
      gids <- region_data[[region_id]]$gid
      res <- compute_region_cor(sids, gids, R_snp, snpinfo$id, weights)
      R_snp <- res$R_snp
      R_snp_gene <- res$R_snp_gene
      R_gene <- res$R_gene
      rownames(R_snp) <- sids
      colnames(R_snp) <- sids
      rownames(R_snp_gene) <- sids
      colnames(R_snp_gene) <- gids
      rownames(R_gene) <- gids
      colnames(R_gene) <- gids
      rm(res)
      # save correlation matrices
      if (isTRUE(save_cor && !is.null(cor_dir))) {
        saveRDS(R_snp_gene, file=R_sg_file)
        saveRDS(R_gene, file=R_g_file)
        saveRDS(R_snp, file=R_s_file)
      }
    }

    finemap_region_res$r2 <- NA
    finemap_region_res$r2[finemap_region_res$type!="SNP"] <- R_gene[finemap_region_res$id[finemap_region_res$type!="SNP"], focus_gid]^2
    finemap_region_res$r2[finemap_region_res$type=="SNP"] <- R_snp_gene[finemap_region_res$id[finemap_region_res$type=="SNP"], focus_gid]^2
    finemap_region_res$r2[finemap_region_res$id == focus_gid] <- 100

    # r2 colors: lead gene: salmon, 0.4~1: purple, others "#7FC97F"
    r2_colors <- c("0-0.4" = "#7FC97F", "0.4-1" = "purple", "1" = "salmon")

    finemap_region_res$r2_levels <- cut(finemap_region_res$r2,
                                        breaks = c(0, 0.4, 1, Inf),
                                        labels = names(r2_colors))
    finemap_region_res$r2_levels <- factor(finemap_region_res$r2_levels, levels = rev(names(r2_colors)))
  } else{
    finemap_region_res$r2 <- NA
    finemap_region_res$r2_levels <- NA
    r2_colors <- "black"
  }

  # gene labels
  finemap_region_res$label <- paste0(finemap_region_res$gene_name, " (", finemap_region_res$context, ")")
  finemap_region_res$label[finemap_region_res$type == "SNP"] <- NA

  # set shapes, sizes, and alpha for data points
  point.shapes <- point.shapes[1:length(unique(finemap_region_res$type))]
  names(point.shapes) <- c("SNP", setdiff(unique(finemap_region_res$type), "SNP"))

  point.alpha <- c("SNP"= point.alpha[1], "non-SNP"= point.alpha[2])
  point.sizes <- c("SNP"= point.sizes[1], "non-SNP"= point.sizes[2])

  # create a locus object for plotting
  loc <- locuszoomr::locus(
    data = finemap_region_res,
    seqname = unique(finemap_region_res$chrom),
    xrange = locus_range,
    ens_db = ens_db,
    labs = "id")

  # get QTL info
  focus_gene_qtls <- rownames(weights[[focus_gid]]$wgt)
  finemap_qtl_res <- finemap_region_res[finemap_region_res$id %in% focus_gene_qtls, ]
  loginfo("QTL positions: %s", finemap_qtl_res$pos)

  # p-value panel
  p_pvalue <- ggplot(loc$data, aes(x = pos/1e6, y = p, label=label,
                                   shape=type, size = object_type, alpha=object_type)) +
    geom_point(aes(color = r2_levels, fill = r2_levels)) +
    geom_text_repel(size=genelabel.size, color="black") +
    scale_shape_manual(values = point.shapes) +
    scale_alpha_manual(values = point.alpha, guide="none") +
    scale_size_manual(values = point.sizes, guide="none") +
    scale_color_manual(values = r2_colors) +
    scale_fill_manual(values = r2_colors, guide="none") +
    xlim(loc$xrange/1e6) +
    labs(x = "", y = expression(-log[10]("p-value")), shape = "Type", color = expression(R^2)) +
    theme_bw() +
    theme(legend.position = legend.position,
          legend.text = element_text(size=legend.text.size),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, l=10, t=10, r=10))

  if (!is.null(p_thresh)) {
    p_pvalue <- p_pvalue + geom_hline(yintercept=-log10(p_thresh), linetype="dashed", color = "red")
  }

  # PIP panel
  p_pip <- ggplot(loc$data, aes(x=pos/1e6, y=susie_pip, label=label,
                                shape=type, size = object_type, alpha=object_type)) +
    geom_point(aes(color = r2_levels, fill = r2_levels)) +
    geom_text_repel(size=genelabel.size, color="black") +
    scale_shape_manual(values = point.shapes) +
    scale_alpha_manual(values = point.alpha, guide="none") +
    scale_size_manual(values = point.sizes, guide="none") +
    scale_color_manual(values = r2_colors) +
    scale_fill_manual(values = r2_colors, guide="none") +
    xlim(loc$xrange/1e6) +
    labs(x = "", y = "PIP", shape = "Type", color = expression(R^2)) +
    theme_bw() +
    theme(legend.position = "none",
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, l=10, t=0, r=10))

  if (!is.null(p_thresh)) {
    p_pip <- p_pip + geom_hline(yintercept=pip_thresh, linetype="dashed", color = "red")
  }

  # QTL panel
  p_qtl <- ggplot(finemap_qtl_res, aes(x=pos/1e6)) +
    geom_rect(aes(xmin=locus_range[1]/1e6, xmax=locus_range[2]/1e6, ymin=0.1, ymax=0.2),
              fill = "lightgray") +
    geom_segment(aes(x=pos/1e6, xend=pos/1e6, y=0.1, yend=0.2), color = "salmon") +
    theme(
      axis.title = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.text.y.left = element_text(angle=0, size=8),
      strip.background = element_blank(),
      plot.margin = margin(b=0, l=10, t=0, r=10))

  # gene track panel
  p_genes <- gg_genetracks(loc, xticks=FALSE, filter_gene_biotype = gene_biotype, text_pos="top")

  if (!is.null(highlight_pos)){
    loginfo("highlight positions: %s", highlight_pos)

    highlight_pos <- as.integer(highlight_pos)

    p_pvalue <- p_pvalue +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dotted", color = "blue", size=0.5)

    p_pip <- p_pip +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dotted", color = "blue", size=0.5)

    p_qtl <- p_qtl +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dotted", color = "blue", size=0.5)

    p_genes <- p_genes +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dotted", color = "blue", size=0.5) +
      theme_bw() +
      theme(legend.position = "none",
            panel.border= element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(b=10, l=10, t=10, r=10))
  }

  cowplot::plot_grid(p_pvalue, p_pip, p_qtl, p_genes, ncol = 1,
                     rel_heights = panel.heights, align = "v")
}
