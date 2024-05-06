
#' Plot the cTWAS result for a single locus
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param region_id region ID to be plotted
#'
#' @param weights a list of weights
#'
#' @param screened_region_data a list of screened region_data for regions with strong signals
#'
#' @param cor_matrix_path path of correlation matrices
#'
#' @param ens_db Ensembl database
#'
#' @param locus_range a vector of start and end positions to define the locus region to be plotted
#'
#' If NULL, the entire region will be plotted (with 100 bp flanks on both sides).
#'
#' @param focus_gene the gene to focus the plot. If NULL, use the lead TWAS feature
#'
#' @param gene_biotype specified biotypes to be displayed in gene tracks or NULL to display all possible ones.
#' By default, only show protein coding genes.
#'
#' @param point.size size values for SNP and non-SNP data points in the scatter plots
#'
#' @param point.alpha alpha values for SNP and non-SNP data points in the scatter plots
#'
#' @param legend.position position to put legends. If "none", no legends will be shown.
#'
#' @param panel.heights Relative heights of the panels.
#'
#' @param check_alignment If TRUE, plots lines at QTL positions in all panels.

#' @importFrom magrittr %>%
#' @import ggplot2
#' @import ggrepel
#' @importFrom locuszoomr locus gg_genetracks
#'
#' @export
#'
make_locusplot <- function(finemap_res,
                           region_id,
                           weights,
                           screened_region_data,
                           cor_matrix_path,
                           ens_db,
                           locus_range = NULL,
                           focus_gene = NULL,
                           gene_biotype = "protein_coding",
                           point.size = c(1, 3.5),
                           point.alpha = c(0.4, 0.6),
                           legend.position = "none",
                           panel.heights = c(4, 4, 0.3, 4),
                           check_alignment = FALSE) {

  # select finemapping result for the target region
  finemap_region_res <- finemap_res[finemap_res$region_id==region_id,]
  finemap_region_res$p <- (-log(2) - pnorm(abs(finemap_region_res$z), lower.tail=F, log.p=T))/log(10) # add pvalue
  finemap_region_res$object_type <- finemap_region_res$type # adding object_type for size and alpha. SNP and non-SNP
  finemap_region_res$object_type[finemap_region_res$object_type!="SNP"] <- "non-SNP"

  if (!is.null(focus_gene)) {
    focus_gene <- finemap_region_res$id[which(finemap_region_res$gene_name == focus_gene)]
  } else{
    # if focus_gene not specified, select the top gene
    finemap_region_gene_res <- finemap_region_res[finemap_region_res$type != "SNP",]
    focus_gene <- finemap_region_gene_res$id[which.max(finemap_region_gene_res$susie_pip)]
  }

  if (!is.null(locus_range)){
    finemap_region_res <- finemap_region_res[finemap_region_res$pos>=locus_range[1] & finemap_region_res$pos<=locus_range[2],, drop=F]
  }else{
    # if locus_range not specified, plot the whole region
    locus_range <- c(min(finemap_region_res$pos)-100, max(finemap_region_res$pos) + 100)
  }

  # add r2 with the focus gene
  R_gene <- readRDS(file.path(cor_matrix_path, paste0("region.", region_id, ".R_gene.RDS")))
  R_snp_gene <- readRDS(file.path(cor_matrix_path, paste0("region.", region_id, ".R_snp_gene.RDS")))
  R_snp <-  readRDS(file.path(cor_matrix_path, paste0("region.", region_id, ".R_snp.RDS")))

  rownames(R_gene) <- screened_region_data[[region_id]]$gid
  colnames(R_gene) <- screened_region_data[[region_id]]$gid
  rownames(R_snp_gene) <- screened_region_data[[region_id]]$sid
  colnames(R_snp_gene) <-screened_region_data[[region_id]]$gid
  rownames(R_snp) <- screened_region_data[[region_id]]$sid
  colnames(R_snp) <- screened_region_data[[region_id]]$sid

  finemap_region_res$r2 <- NA
  finemap_region_res$r2[finemap_region_res$type!="SNP"] <- R_gene[finemap_region_res$id[finemap_region_res$type!="SNP"], focus_gene]^2
  finemap_region_res$r2[finemap_region_res$type=="SNP"] <- R_snp_gene[finemap_region_res$id[finemap_region_res$type=="SNP"], focus_gene]^2
  finemap_region_res$r2[finemap_region_res$id == focus_gene] <- 100

  # r2 colors: lead gene: salmon, 0.4~1: purple, others "#7FC97F"
  r2_colors <- c("0-0.4" = "#7FC97F", "0.4-1" = "purple", "1" = "salmon")

  finemap_region_res$r2_levels <- cut(finemap_region_res$r2,
                                      breaks = c(0, 0.4, 1, Inf),
                                      labels = names(r2_colors))
  finemap_region_res$r2_levels <- factor(finemap_region_res$r2_levels, levels = rev(names(r2_colors)))

  # gene labels
  finemap_region_res$label <- paste0(finemap_region_res$gene_name, "(", finemap_region_res$context, ")")
  finemap_region_res$label[finemap_region_res$type == "SNP"] <- NA

  # create a locus object for plotting
  loc <- locuszoomr::locus(
    data = finemap_region_res,
    seqname = unique(finemap_region_res$chrom),
    xrange = locus_range,
    ens_db = ens_db,
    labs = "id")

  # loc$data$r2_color <- finemap_region_res$r2_color

  # get QTL info
  focus_gene_qtls <- rownames(weights[[focus_gene]]$wgt)
  finemap_qtl_res <- finemap_region_res[finemap_region_res$id %in% focus_gene_qtls, ]
  # focus_gene_QTLs_res$qtl_name <- NULL

  # p-value panel
  p_pvalue <- ggplot(loc$data, aes(x = pos/1e6, y = p, label=label)) +
    geom_point(aes(shape=type, size = object_type, alpha=object_type, color = r2_levels, fill = r2_levels)) +
    geom_text_repel(size=2, color="black") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    scale_alpha_manual(values=c("SNP"=point.alpha[1], "non-SNP"=point.alpha[2]), guide="none") +
    scale_size_manual(values=c("SNP"=point.size[1], "non-SNP"=point.size[2]), guide="none") +
    scale_color_manual(values = r2_colors) +
    scale_fill_manual(values = r2_colors, guide="none") +
    xlim(loc$xrange/1e6) +
    labs(x = "", y = expression(-log[10]("p-value")), shape = "Type", color = expression(R^2)) +
    theme_bw() +
    theme(legend.position = legend.position,
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, l=10, t=10, r=10))


  # PIP panel
  p_pip <- ggplot(loc$data, aes(x=pos/1e6, y=susie_pip, label=label)) +
    geom_point(aes(shape=type, size = object_type, alpha=object_type, color = r2_levels, fill = r2_levels)) +
    geom_text_repel(size=2, color="black") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    scale_alpha_manual(values=c("SNP"=point.alpha[1], "non-SNP"=point.alpha[2]), guide="none") +
    scale_size_manual(values=c("SNP"=point.size[1], "non-SNP"=point.size[2]), guide="none") +
    scale_color_manual(values = r2_colors) +
    scale_fill_manual(values = r2_colors, guide="none") +
    xlim(loc$xrange/1e6) +
    labs(x = "", y = "PIP", shape = "Type", color = expression(R^2)) +
    theme_bw() +
    theme(legend.position = "none",
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, l=10, t=0, r=10))

  # QTL panel
  p_qtl <- ggplot(finemap_qtl_res) +
    geom_rect(aes(xmin=locus_range[1]/1e6, xmax=locus_range[2]/1e6, ymin=0.1, ymax=0.2),
              fill = "lightgray") +
    geom_segment(aes(x=pos/1e6, xend=pos/1e6, y=0.1, yend=0.2), color = r2_colors["1"]) +
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

  if (check_alignment){
    p_pvalue <- p_pvalue +
      geom_vline(xintercept = finemap_qtl_res$pos/1e6, linetype="dotted", color = "blue", size=0.5)

    p_pip <- p_pip +
      geom_vline(xintercept = finemap_qtl_res$pos/1e6, linetype="dotted", color = "blue", size=0.5)

    p_qtl <- p_qtl +
      geom_vline(xintercept = finemap_qtl_res$pos/1e6, linetype="dotted", color = "blue", size=0.5)

    p_genes <- p_genes +
      geom_vline(xintercept = finemap_qtl_res$pos/1e6, linetype="dotted", color = "blue", size=0.5) +
      theme_bw() +
      theme(legend.position = "none",
            panel.border= element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(b=10, l=10, t=10, r=10))
  }

  cowplot::plot_grid(p_pvalue, p_pip, p_qtl, p_genes, ncol = 1,
                     rel_heights = panel.heights, align = "v")
}
