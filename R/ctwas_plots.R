
#' @title Plots the cTWAS result for a single locus
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param region_id region ID to be plotted
#'
#' @param weights a list of weights
#'
#' @param ens_db Ensembl database
#'
#' @param R_snp_gene SNP-gene correlation matrix of the region.
#' If both R_snp_gene and R_gene are available,
#' color data points with correlations with the focus gene.
#'
#' @param R_gene gene-gene correlation matrix of the region
#'
#' @param locus_range a vector of start and end positions to define the region boundary to be plotted.
#' If NULL, the entire region will be plotted (with 100 bp flanks on both sides).
#'
#' @param focus_gene the gene name to focus the plot.
#' If NULL, choose the gene with the highest PIP.
#'
#' @param filter_gene_biotype biotype to be displayed in gene tracks.
#' By default, limits to protein coding genes. If NULL, display all possible ones.
#'
#' @param highlight_pval pvalue to highlight with a horizontal line
#'
#' @param highlight_pip PIP to highlight with a horizontal line
#'
#' @param highlight_pos genomic positions to highlight with vertical lines
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
#' @importFrom magrittr %>%
#' @importFrom locuszoomr locus gg_genetracks
#' @importFrom logging loginfo
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 margin
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
make_locusplot <- function(finemap_res,
                           region_id,
                           weights,
                           ens_db,
                           R_snp_gene = NULL,
                           R_gene = NULL,
                           locus_range = NULL,
                           focus_gene = NULL,
                           filter_gene_biotype = "protein_coding",
                           highlight_pval = NULL,
                           highlight_pip = 0.8,
                           highlight_pos = NULL,
                           point.sizes = c(1, 3.5),
                           point.alpha = c(0.4, 0.6),
                           point.shapes = c(21:25),
                           genelabel.size = 2.5,
                           legend.text.size = 10,
                           legend.position = c("none", "top", "left", "right", "bottom"),
                           panel.heights = c(4, 4, 0.3, 4)) {

  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }

  legend.position <- match.arg(legend.position)

  # select finemapping result for the target region
  finemap_region_res <- finemap_res[which(finemap_res$region_id==region_id), ]
  # convert z to -log10(pval)
  finemap_region_res$p <- (-log(2) - pnorm(abs(finemap_region_res$z), lower.tail=F, log.p=T))/log(10)
  # add object_type to plot SNP and non-SNP categories with different sizes and alphas
  finemap_region_res$object_type <- finemap_region_res$type
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
    # if focus_gene is not specified, choose the top gene with highest PIP
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
    finemap_region_res <- finemap_region_res[idx_in_range,, drop=FALSE]
  }else{
    # if locus_range not specified, plot the whole region
    locus_range <- c(min(finemap_region_res$pos)-100, max(finemap_region_res$pos) + 100)
  }
  loginfo("locus_range: %s", locus_range)

  if (!is.null(R_gene) && !is.null(R_snp_gene)) {
    plot_r2 <- TRUE
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
    plot_r2 <- FALSE
    # finemap_region_res$r2 <- NA
    # finemap_region_res$r2_levels <- NA
    # r2_colors <- "black"
    r2_colors <- c("0-1" = "gray", "1" = "salmon")

    finemap_region_res$r2 <- 0.1
    finemap_region_res$r2[finemap_region_res$id == focus_gid] <- 100
    finemap_region_res$r2_levels <- cut(finemap_region_res$r2,
                                        breaks = c(0, 1, Inf),
                                        labels = names(r2_colors))
    finemap_region_res$r2_levels <- factor(finemap_region_res$r2_levels, levels = rev(names(r2_colors)))
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
  loc <- locus(
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

  if (plot_r2) {
    p_pvalue <- p_pvalue +
      scale_color_manual(values = r2_colors) +
      scale_fill_manual(values = r2_colors, guide="none")
  } else {
    p_pvalue <- p_pvalue +
      scale_color_manual(values = r2_colors, guide="none") +
      scale_fill_manual(values = r2_colors, guide="none")
  }

  if (!is.null(highlight_pval)) {
    p_pvalue <- p_pvalue +
      geom_hline(yintercept=-log10(highlight_pval), linetype="dashed", color = "red")
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

  if (!is.null(highlight_pip)) {
    p_pip <- p_pip +
      geom_hline(yintercept=highlight_pip, linetype="dashed", color = "red")
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
  p_genes <- gg_genetracks(loc, xticks=FALSE,
                           filter_gene_biotype = filter_gene_biotype,
                           text_pos="top")

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
                     rel_heights = panel.heights, align = "v", axis="tblr")
}
