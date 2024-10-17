
#' @title Plots the cTWAS result for a single locus
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#'
#' @param region_id region ID to be plotted
#'
#' @param ens_db Ensembl database
#'
#' @param weights a list of proprocessed weights, used to plot QTL track.
#' If NULL, will not plot the QTL track.
#'
#' @param R_snp_gene SNP-gene correlation matrix of the region.
#' If both R_snp_gene and R_gene are available,
#' color data points with correlations with the focal gene.
#'
#' @param R_gene gene-gene correlation matrix of the region
#'
#' @param locus_range a vector of start and end positions to define the region boundary to be plotted.
#' If NULL, the entire region will be plotted (with 100 bp flanks on both sides).
#'
#' @param focal_id the focal ID.
#'
#' @param focal_gene the focal gene name. By default, use the gene with the highest PIP.
#'
#' @param filter_protein_coding_genes If TRUE, limits to protein coding genes only.
#'
#' @param filter_cs If TRUE, limits to credible sets.
#'
#' @param color_pval_by Options to color the p-value track.
#' "LD": colors the p-value track by the correlations to the focal gene.
#' "cs": colors the p-value track by the credible sets.
#' "none": uses the same color for non-focal genes.
#'
#' @param color_pip_by Options to color the PIP track.
#' "LD": colors the PIP track by the correlations to the focal gene.
#' "cs": colors the PIP track by the credible sets.
#' "none": uses the same color for non-focal genes.
#'
#' @param LD.colors Colors for correlation levels.
#'
#' @param cs.colors Colors for credible sets.
#'
#' @param focal.colors Colors for non-focal and focal gene.
#'
#' @param label_QTLs If TRUE, label SNP IDs in the QTL panel.
#'
#' @param highlight_pval p-value to highlight with a horizontal line
#'
#' @param highlight_pip PIP to highlight with a horizontal line
#'
#' @param highlight_pos genomic positions to highlight with vertical lines
#'
#' @param highlight.color color for horizontal lines highlighting p-values and PIPs.
#'
#' @param highlight.pos.color color for vertical lines highlighting genomic positions.
#'
#' @param point.sizes size values for SNP and non-SNP data points in the scatter plots
#'
#' @param point.alpha alpha values for SNP and non-SNP data points in the scatter plots
#'
#' @param point.shapes shapes values for data points of different types in the scatter plots
#'
#' @param label.text.size Font size for gene and SNP label text
#'
#' @param max.overlaps Setting for geom_text_repel() function to label texts.
#' Exclude text labels when they overlap too many other things.
#'
#' @param axis.text.size Font size for axis label text.
#'
#' @param axis.title.size Font size for axis title text.
#'
#' @param legend.text.size Font size for legend text.
#'
#' @param legend.position position to put legends. If "none", no legends will be shown.
#'
#' @param genelabel.cex.text Size for gene label text.
#'
#' @param panel.heights Relative heights of the panels.
#'
#' @param verbose If TRUE, print detail messages.
#'
#' @importFrom magrittr %>%
#' @importFrom locuszoomr locus gg_genetracks
#' @importFrom logging loginfo logwarn
#' @importFrom cowplot plot_grid theme_cowplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_brewer
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
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
#' @importFrom rlang .data
#' @importFrom grid unit
#'
#' @export
#'
make_locusplot <- function(finemap_res,
                           region_id,
                           ens_db,
                           weights = NULL,
                           R_snp_gene = NULL,
                           R_gene = NULL,
                           locus_range = NULL,
                           focal_id = NULL,
                           focal_gene = NULL,
                           filter_protein_coding_genes = TRUE,
                           filter_cs = TRUE,
                           color_pval_by = c("cs", "LD", "none"),
                           color_pip_by = c("cs", "LD", "none"),
                           LD.colors = c("grey", "blue", "purple", "salmon"),
                           cs.colors = c("firebrick", "dodgerblue", "forestgreen", "darkmagenta", "darkorange"),
                           focal.colors = c("grey", "salmon"),
                           label_QTLs = TRUE,
                           highlight_pval = NULL,
                           highlight_pip = 0.8,
                           highlight_pos = NULL,
                           highlight.color = "red",
                           highlight.pos.color = "cyan2",
                           point.sizes = c(1, 3),
                           point.alpha = c(0.4, 0.6),
                           point.shapes = c(16, 15, 18, 17, 10, 12, 14, 11),
                           label.text.size = 2.5,
                           max.overlaps = 10,
                           axis.text.size = 8,
                           axis.title.size = 10,
                           legend.text.size = 9,
                           legend.position = "top",
                           genelabel.cex.text = 0.7,
                           panel.heights = c(4, 4, 1, 4),
                           verbose = FALSE) {

  color_pval_by <- match.arg(color_pval_by)
  color_pip_by <- match.arg(color_pip_by)

  if (!inherits(ens_db,"EnsDb"))
    stop("'ens_db' should be a EnsDb object.")

  if (is.null(finemap_res$chrom) || is.null(finemap_res$pos))
    stop("Please add 'chrom' and 'pos' columns to finemapping result!")

  if (anyNA(finemap_res$pos))
    stop("Missing values found in the 'pos' column of finemapping result!")

  # input data should be a data frame
  finemap_res <- as.data.frame(finemap_res)
  # select finemapping result for the target region
  finemap_region_res <- finemap_res[which(finemap_res$region_id==region_id), ]

  if (!"pval" %in% names(finemap_region_res)){
    # convert z to -log10(pval)
    finemap_region_res$pval <- z2p(finemap_region_res$z, neg_log10_p = TRUE)
  } else {
    if (max(finemap_region_res$pval) <= 1) {
      finemap_region_res$pval <- -log10(finemap_region_res$pval)
    }
  }

  # add object_type to plot SNP and non-SNP categories with different sizes and alphas
  finemap_region_res$object_type <- finemap_region_res$type
  finemap_region_res$object_type[finemap_region_res$object_type!="SNP"] <- "non-SNP"

  # gene labels
  # if gene_name is not in finemap_res, use molecular_id.
  if (is.null(finemap_region_res$gene_name)){
    loginfo("'gene_name' not found in finemap_res. Use 'molecular_id' instead.")
    if (is.null(finemap_region_res$molecular_id)) {
      finemap_region_res$molecular_id <- get_molecular_ids(finemap_region_res)
    }
    finemap_region_res$gene_name <- finemap_region_res$molecular_id
  }
  finemap_region_res$label <- paste0(finemap_region_res$gene_name, " (", finemap_region_res$context, ")")
  finemap_region_res$label[finemap_region_res$group == "SNP"] <- NA

  # limit results and gene tracks to protein coding genes
  if (filter_protein_coding_genes) {
    if (is.null(finemap_region_res$gene_type)) {
      loginfo("'gene_type' column cannot be found in finemap_res. Skipped filtering protein coding genes.")
    } else {
      loginfo("Limit to protein coding genes")
      drop_idx <- which(finemap_region_res$group!="SNP" & finemap_region_res$gene_type!="protein_coding")
      if (length(drop_idx) > 0) {
        finemap_region_res <- finemap_region_res[-drop_idx, , drop=FALSE]
      }
    }
    filter_gene_biotype <- "protein_coding"
  } else{
    filter_gene_biotype <- NULL
  }

  # If focal_gene is specified, extract the focal id;
  # otherwise, choose the top gene with highest PIP
  if (is.null(focal_id)) {
    if (!is.null(focal_gene)) {
      focal_idx <- which(finemap_region_res$gene_name == focal_gene)
      if (length(focal_idx) == 0){
        stop("can't find focal gene in the fine-mapping result of the region")
      }
      if (length(focal_idx) > 1){
        finemap_region_focal_gene_res <- finemap_region_res[focal_idx,]
        focal_id <- finemap_region_focal_gene_res$id[which.max(finemap_region_focal_gene_res$susie_pip)]
      } else {
        focal_id <- finemap_region_res$id[focal_idx]
      }
    } else {
      # if focal_gene is not specified, choose the top gene with highest PIP
      finemap_region_gene_res <- finemap_region_res[finemap_region_res$group != "SNP",]
      focal_idx <- which.max(finemap_region_gene_res$susie_pip)
      focal_gene <- finemap_region_gene_res$gene_name[focal_idx]
      focal_id <- finemap_region_gene_res$id[focal_idx]
    }
  }

  loginfo("focal id: %s", focal_id)
  focal_gene_name <- finemap_region_res$gene_name[finemap_region_res$id == focal_id]
  focal_gene_context <- finemap_region_res$context[finemap_region_res$id == focal_id]
  focal_gene_type <- finemap_region_res$type[finemap_region_res$id == focal_id]
  loginfo("focal molecular trait: %s %s %s", focal_gene_name, focal_gene_context, focal_gene_type)

  # set locus range for plotting
  if (!is.null(locus_range)){
    idx_in_range <- which(finemap_region_res$pos>=locus_range[1] & finemap_region_res$pos<=locus_range[2])
    finemap_region_res <- finemap_region_res[idx_in_range,, drop=FALSE]
  }else{
    # if locus_range not specified, plot the whole region
    locus_range <- c(min(finemap_region_res$pos)-100, max(finemap_region_res$pos) + 100)
  }
  chrom <- unique(finemap_region_res$chrom)
  loginfo("Range of locus: chr%s:%s-%s", chrom, locus_range[1], locus_range[2])

  # set colors for LD with focal gene
  if (color_pval_by == "LD" || color_pip_by == "LD") {
    if (is.null(R_gene) || is.null(R_snp_gene)) {
      loginfo("'R_gene' or 'R_snp_gene' not available. Cannot color by LD!")
      color_pval_by[color_pval_by == "LD"] <- "none"
      color_pip_by[color_pip_by == "LD"] <- "none"
    } else {
      finemap_region_res$r2 <- NA
      finemap_region_res$r2[finemap_region_res$type!="SNP"] <- R_gene[finemap_region_res$id[finemap_region_res$type!="SNP"], focal_id]^2
      finemap_region_res$r2[finemap_region_res$type=="SNP"] <- R_snp_gene[finemap_region_res$id[finemap_region_res$type=="SNP"], focal_id]^2
      finemap_region_res$r2[finemap_region_res$id == focal_id] <- 100
      LD_colors <- c("0-0.2" = LD.colors[1], "0.2-0.4" = LD.colors[2],
                     "0.4-1" = LD.colors[3], "1" = LD.colors[4])
      finemap_region_res$r2_levels <- cut(finemap_region_res$r2,
                                          breaks = c(0, 0.2, 0.4, 1, Inf),
                                          labels = names(LD_colors))
      finemap_region_res$r2_levels <- factor(finemap_region_res$r2_levels, levels = rev(names(LD_colors)))
    }
  }

  # set colors for credible sets
  if (color_pval_by == "cs" || color_pip_by == "cs") {
    if (is.null(finemap_region_res$cs)){
      logwarn("'cs' not available. Cannot coloring by cs!")
      color_pval_by[color_pval_by == "cs"] <- "none"
      color_pip_by[color_pip_by == "cs"] <- "none"
    } else {
      cs_colors <- c("L1" = cs.colors[1], "L2" = cs.colors[2], "L3" = cs.colors[3], "L4" = cs.colors[4], "L5" = cs.colors[5])
      # if there are multiple CSs, use the first one
      if (any(grepl(",", finemap_region_res$cs))){
        finemap_region_res$cs <- sapply(strsplit(finemap_region_res$cs, ","), "[[", 1)
      }
      finemap_region_res$cs <- factor(finemap_region_res$cs, levels = names(cs_colors))
    }
  }

  # set colors for focal and non-focal genes
  if (color_pval_by == "none" || color_pip_by == "none") {
    focal_colors <- c("non-focal" = focal.colors[1], "focal" = focal.colors[2])
    finemap_region_res$focal_levels <- "non-focal"
    finemap_region_res$focal_levels[finemap_region_res$id == focal_id] <- "focal"
  }

  # set shapes, sizes, and alpha for data points
  finemap_region_res$group <- factor(finemap_region_res$group,
                                     levels = c(setdiff(unique(finemap_region_res$group), "SNP"), "SNP"))

  finemap_region_res$type <- factor(finemap_region_res$type,
                                    levels = c(setdiff(unique(finemap_region_res$type), "SNP"), "SNP"))

  finemap_region_res$context <- factor(finemap_region_res$context,
                                       levels = c(setdiff(unique(finemap_region_res$context), "SNP"), "SNP"))

  finemap_region_res$object_type <- factor(finemap_region_res$object_type,
                                           levels = c("non-SNP", "SNP"))

  region_types <- levels(finemap_region_res$type)
  region_contexts <- levels(finemap_region_res$context)
  region_group <- levels(finemap_region_res$group)
  region_nonSNP_types <- setdiff(region_types, "SNP")

  point.shapes <- c(point.shapes[2:length(region_types)], point.shapes[1])
  names(point.shapes) <- c(region_nonSNP_types, "SNP")

  point.alpha <- c("non-SNP" = point.alpha[2], "SNP" = point.alpha[1])
  point.sizes <- c("non-SNP" = point.sizes[2], "SNP" = point.sizes[1])

  legend.sizes <- c(rep(3.5, length(region_nonSNP_types)), 2)
  names(legend.sizes) <- c(region_nonSNP_types, "SNP")

  # create a locus object for plotting
  loc <- locus(
    data = finemap_region_res,
    seqname = chrom,
    xrange = locus_range,
    ens_db = ens_db,
    labs = "id")

  # get QTLs for the focal gene
  if (!is.null(weights)){
    if (!inherits(weights,"list"))
      stop("'weights' should be a list!")

    focal_gene_qtls <- rownames(weights[[focal_id]]$wgt)
    finemap_qtl_res <- finemap_region_res[finemap_region_res$id %in% focal_gene_qtls, ]
    loginfo("focal molecular trait QTL positions: %s", finemap_qtl_res$pos)
  }

  # p-value panel
  if (verbose) {
    loginfo("Making p-value panel ...")
  }
  pval_plot_data <- loc$data
  p_pval <- ggplot(pval_plot_data, aes(x=.data$pos/1e6, y=.data$pval, shape=.data$type,
                                       size=.data$object_type, alpha=.data$object_type)) +
    geom_text_repel(aes(label=.data$label), size=label.text.size, color="black",
                    max.overlaps = max.overlaps, na.rm = TRUE) +
    scale_shape_manual(values = point.shapes) +
    scale_alpha_manual(values = point.alpha, guide="none") +
    scale_size_manual(values = point.sizes, guide="none") +
    xlim(loc$xrange/1e6) +
    labs(x = "", y = expression(-log[10]("p-val")), shape = "") +
    theme_bw() +
    theme(legend.position = legend.position,
          legend.spacing.x = unit(1.0, 'cm'),
          legend.title = element_text(size=legend.text.size),
          legend.text = element_text(size=legend.text.size),
          axis.title = element_text(size=axis.title.size),
          axis.text = element_text(size=axis.text.size),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, l=10, t=10, r=10))

  if (color_pval_by == "cs") {
    p_pval <- p_pval +
      geom_point(aes(color=.data$cs)) +
      scale_color_manual(values = cs_colors) +
      labs(color = "CS") +
      guides(shape = guide_legend(order = 1, override.aes = list(size = legend.sizes)),
             color = guide_legend(order = 2, override.aes = list(size = 1.5)))
  } else if (color_pval_by == "LD") {
    p_pval <- p_pval +
      geom_point(aes(color=.data$r2_levels)) +
      scale_color_manual(values = LD_colors) +
      labs(color = expression(R^2)) +
      guides(shape = guide_legend(order = 1, override.aes = list(size = legend.sizes)),
             color = guide_legend(order = 2, override.aes = list(size = 1.5)))
  } else {
    p_pval <- p_pval +
      geom_point(aes(color=.data$focal_levels)) +
      scale_color_manual(values = focal_colors, guide="none") +
      guides(shape = guide_legend(override.aes = list(size = legend.sizes)))
  }

  if (!is.null(highlight_pval)) {
    p_pval <- p_pval +
      geom_hline(yintercept=-log10(highlight_pval), linetype="dashed", color = highlight.color)
  }

  # PIP panel
  if (verbose) {
    loginfo("Making PIP panel ...")
  }
  pip_plot_data <- loc$data
  # limit to credible sets (if cs is available)
  if (filter_cs && !is.null(pip_plot_data$cs)) {
    loginfo("Limit PIPs to credible sets")
    pip_plot_data <- pip_plot_data[!is.na(pip_plot_data$cs),]
  }

  p_pip <- ggplot(pip_plot_data, aes(x=.data$pos/1e6, y=.data$susie_pip, shape=.data$type,
                                     size=.data$object_type, alpha=.data$object_type)) +
    geom_text_repel(aes(label=.data$label), size=label.text.size, color="black", max.overlaps = max.overlaps, na.rm = TRUE) +
    scale_shape_manual(values = point.shapes, guide="none") +
    scale_alpha_manual(values = point.alpha, guide="none") +
    scale_size_manual(values = point.sizes, guide="none") +
    xlim(loc$xrange/1e6) +
    ylim(0,1) +
    labs(x = "", y = "cTWAS PIP", shape = "") +
    theme_bw() +
    theme(legend.position = legend.position,
          legend.spacing.x = unit(1.0, 'cm'),
          legend.title = element_text(size=legend.text.size),
          legend.text = element_text(size=legend.text.size),
          axis.title = element_text(size=axis.title.size),
          axis.text = element_text(size=axis.text.size),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, l=10, t=10, r=10))

  if (color_pip_by == "cs") {
    p_pip <- p_pip +
      geom_point(aes(color=.data$cs)) +
      scale_color_manual(values = cs_colors) +
      labs(color = "CS") +
      guides(color = guide_legend(override.aes = list(size = 1.5)))
  } else if (color_pip_by == "LD") {
    p_pip <- p_pip +
      geom_point(aes(color=.data$r2_levels)) +
      scale_color_manual(values = LD_colors) +
      labs(color = expression(R^2)) +
      guides(color = guide_legend(override.aes = list(size = 1.5)))
  } else {
    p_pip <- p_pip +
      geom_point(aes(color=.data$focal_levels)) +
      scale_color_manual(values = focal_colors, guide="none")
  }

  if (!is.null(highlight_pip)) {
    p_pip <- p_pip +
      geom_hline(yintercept=highlight_pip, linetype="dashed", color = highlight.color)
  }

  # QTL panel
  if (!is.null(weights)){
    if (verbose) {
      loginfo("Making QTL panel ...")
    }
    p_qtl <- ggplot(finemap_qtl_res, aes(x=.data$pos/1e6)) +
      geom_rect(aes(xmin=loc$xrange[1]/1e6, xmax=loc$xrange[2]/1e6, ymin=0, ymax=1),
                fill="gray90") +
      geom_segment(aes(x=.data$pos/1e6, xend=.data$pos/1e6, y=0, yend=1), color=focal.colors[2]) +
      labs(title = paste(focal_gene_name, focal_gene_context, focal_gene_type),
           y = "QTL") +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_text(size=axis.title.size),
            axis.text = element_text(size=axis.text.size),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size=legend.text.size),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(b=0, l=10, t=0, r=10))

    if (label_QTLs){
      p_qtl <- p_qtl +
        geom_text_repel(aes(y=0.5, label=.data$id), size=label.text.size,
                        color="black", max.overlaps = max.overlaps)
    }
  } else {
    p_qtl <- NULL
  }

  # gene track panel
  if (verbose) {
    loginfo("Making gene track panel ...")
  }
  p_genes <- gg_genetracks(loc,
                           filter_gene_biotype = filter_gene_biotype,
                           xticks = FALSE,
                           text_pos = "top",
                           cex.text = genelabel.cex.text) +
    labs(x = paste0("chr", chrom)) +
    theme_bw() +
    theme(legend.position = "none",
          panel.border = element_blank(),
          axis.title = element_text(size=axis.text.size+1),
          axis.text = element_text(size=axis.text.size),
          plot.margin = margin(b=10, l=10, t=0, r=10))

  if (!is.null(highlight_pos)){
    loginfo("highlight positions: %s", highlight_pos)
    highlight_pos <- as.integer(highlight_pos)

    p_pval <- p_pval +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dashed", color = highlight.pos.color)

    p_pip <- p_pip +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dashed", color = highlight.pos.color)

    p_qtl <- p_qtl +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dashed", color = highlight.pos.color)

    p_genes <- p_genes +
      geom_vline(xintercept = highlight_pos/1e6, linetype="dashed", color = highlight.pos.color)
  }

  if (verbose) {
    loginfo("Plot panels ...")
  }

  if (is.null(p_qtl)){
    panel.heights <- panel.heights[c(1,2,length(panel.heights))]
    cowplot::plot_grid(p_pval, p_pip, p_genes, ncol = 1,
                       rel_heights = panel.heights, align = "v", axis="tblr")
  } else {
    cowplot::plot_grid(p_pval, p_pip, p_qtl, p_genes, ncol = 1,
                       rel_heights = panel.heights, align = "v", axis="tblr")
  }

}



#' @title Make convergence plots for the estimated parameters
#'
#' @param param a list of cTWAS parameter estimation result from \code{est_param}
#'
#' @param gwas_n the sample size of the GWAS summary statistics
#'
#' @param title.size font size of the plot title
#'
#' @param legend.size font size of the plot legend title
#'
#' @param colors colors for different groups, passed as the
#'   \code{values} input to \code{\link[ggplot2]{scale_color_manual}}.
#'   If fewer colors than "fits" are given, the colors are recycled.
#'
#' @param ncol number of columns in the output plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab ylab ggtitle
#' @importFrom ggplot2 expand_limits
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 scale_color_manual
#' @importFrom cowplot plot_grid theme_cowplot
#' @importFrom rlang .data
#'
#' @export
make_convergence_plots <- function(param,
                                   gwas_n,
                                   title.size = 10,
                                   legend.size = 8,
                                   colors = c("#E69F00","#56B4E9","#009E73","#F0E442",
                                              "#0072B2","#D55E00","#CC79A7", "#999999"),
                                   ncol = 2){

  # estimated group prior (all iterations)
  group_prior_iters <- param$group_prior_iters
  # estimated group prior variance (all iterations)
  group_prior_var_iters <- param$group_prior_var_iters

  # group size
  group_size <- param$group_size[rownames(group_prior_iters)]
  group_size <- as.numeric(group_size)
  names(group_size) <- rownames(group_prior_iters)

  # estimated group PVE (all iterations)
  group_pve_iters <- group_prior_var_iters*group_prior_iters*group_size/gwas_n

  # estimated enrichment of genes (all iterations)
  enrichment_iters <- t(sapply(rownames(group_prior_iters)[rownames(group_prior_iters)!="SNP"], function(x){
    group_prior_iters[rownames(group_prior_iters)==x,]/group_prior_iters[rownames(group_prior_iters)=="SNP"]}))

  # make convergence plots

  # prior inclusion probability plot
  df <- data.frame(niter = rep(1:ncol(group_prior_iters), nrow(group_prior_iters)),
                   value = unlist(lapply(1:nrow(group_prior_iters), function(x){group_prior_iters[x,]})),
                   group = rep(rownames(group_prior_iters), each=ncol(group_prior_iters)))
  factor_levels <- c(setdiff(rownames(group_prior_iters), "SNP"), "SNP")
  df$group <- factor(df$group, levels = factor_levels)
  df <- na.omit(df)

  p_pi <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab(bquote(pi)) +
    ggtitle("Proportion Causal") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size),
          plot.margin = margin(b=10, l=10, t=10, r=10))

  # prior effect size variance plot
  df <- data.frame(niter = rep(1:ncol(group_prior_var_iters), nrow(group_prior_var_iters)),
                   value = unlist(lapply(1:nrow(group_prior_var_iters), function(x){group_prior_var_iters[x,]})),
                   group = rep(rownames(group_prior_var_iters), each=ncol(group_prior_var_iters)))
  df$group <- factor(df$group, levels = factor_levels)
  df <- na.omit(df)

  p_sigma2 <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab(bquote(sigma^2)) +
    ggtitle("Effect Size") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size),
          plot.margin = margin(b=10, l=10, t=10, r=10))

  # enrichment plot
  df <- data.frame(niter = rep(1:ncol(enrichment_iters), nrow(enrichment_iters)),
                   value = unlist(lapply(1:nrow(enrichment_iters), function(x){enrichment_iters[x,]})),
                   group = rep(rownames(enrichment_iters), each=ncol(enrichment_iters)))
  df$group <- factor(df$group, levels = factor_levels[factor_levels!="SNP"])
  df <- na.omit(df)

  p_enrich <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab(bquote(pi[G]/pi[S])) +
    ggtitle("Enrichment") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size),
          plot.margin = margin(b=10, l=10, t=10, r=10))

  # PVE plot
  df <- data.frame(niter = rep(1:ncol(group_pve_iters), nrow(group_pve_iters)),
                   value = unlist(lapply(1:nrow(group_pve_iters), function(x){group_pve_iters[x,]})),
                   group = rep(rownames(group_pve_iters), each=ncol(group_pve_iters)))
  df$group <- factor(df$group, levels = factor_levels)
  df <- na.omit(df)

  p_pve <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab("PVE") +
    ggtitle("PVE") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size),
          plot.margin = margin(b=10, l=10, t=10, r=10))

  cowplot::plot_grid(p_pi, p_sigma2, p_enrich, p_pve,
                     ncol = ncol)
}
