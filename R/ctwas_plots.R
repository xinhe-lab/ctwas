
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
#' @param filter_genetrack_biotype biotype to be displayed in gene tracks.
#' By default, limits to protein coding genes. If NULL, display all possible ones.
#'
#' @param highlight_pval pvalue to highlight with a horizontal line
#'
#' @param highlight_pip PIP to highlight with a horizontal line
#'
#' @param highlight_pos genomic positions to highlight with vertical lines
#'
#' @param label_QTLs If TRUE, label SNP IDs in the QTL panel.
#'
#' @param point.sizes size values for SNP and non-SNP data points in the scatter plots
#'
#' @param point.alpha alpha values for SNP and non-SNP data points in the scatter plots
#'
#' @param point.shapes shapes values for data points of different types in the scatter plots
#'
#' @param label.text.size Font size for gene and SNP label text
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
#' @importFrom cowplot plot_grid theme_cowplot
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
#' @importFrom rlang .data
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
                           filter_genetrack_biotype = "protein_coding",
                           highlight_pval = NULL,
                           highlight_pip = 0.8,
                           highlight_pos = NULL,
                           label_QTLs = TRUE,
                           point.sizes = c(1, 3.5),
                           point.alpha = c(0.4, 0.6),
                           point.shapes = c(21:25),
                           label.text.size = 2.5,
                           legend.text.size = 10,
                           legend.position = "top",
                           panel.heights = c(4, 4, 1, 4)) {

  if (!inherits(weights,"list")){
    stop("'weights' should be a list.")
  }

  # Check to see if gene_name and gene_type are already in finemap_res
  annot_cols <- c("gene_name", "gene_type")
  if (!all(annot_cols %in% colnames(finemap_res))){
    stop("finemap_res needs to contain the following columns: ",
         paste(annot_cols, collapse = " "),
         "\nPlease first run anno_finemap_res() to annotate finemap_res")
  }

  # input data should be a data frame
  finemap_res <- as.data.frame(finemap_res)
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
  chrom <- unique(finemap_region_res$chrom)
  loginfo("plot locus range: chr%s %s", chrom, locus_range)

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
    r2_colors <- c("0-1" = "gray50", "1" = "salmon")
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
    seqname = chrom,
    xrange = locus_range,
    ens_db = ens_db,
    labs = "id")

  # get QTLs for the focus gene
  focus_gene_qtls <- rownames(weights[[focus_gid]]$wgt)
  finemap_qtl_res <- finemap_region_res[finemap_region_res$id %in% focus_gene_qtls, ]
  focus_gene_name <- finemap_region_res$gene_name[finemap_region_res$id == focus_gid]
  focus_gene_context <- finemap_region_res$context[finemap_region_res$id == focus_gid]
  focus_gene_type <- finemap_region_res$type[finemap_region_res$id == focus_gid]
  loginfo("%s %s %s QTLs", focus_gene_name, focus_gene_context, focus_gene_type)
  loginfo("QTL positions: %s", finemap_qtl_res$pos)

  # p-value panel
  p_pvalue <- ggplot(loc$data, aes(x=.data$pos/1e6, y=.data$p, shape=.data$type,
                                   size=.data$object_type, alpha=.data$object_type)) +
    geom_point(aes(color=.data$r2_levels, fill=.data$r2_levels)) +
    geom_text_repel(aes(label=.data$label), size=label.text.size, color="black") +
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
  p_pip <- ggplot(loc$data, aes(x=.data$pos/1e6, y=.data$susie_pip, shape=.data$type,
                                size=.data$object_type, alpha=.data$object_type)) +
    geom_point(aes(color=.data$r2_levels, fill=.data$r2_levels)) +
    geom_text_repel(aes(label=.data$label), size=label.text.size, color="black") +
    scale_shape_manual(values = point.shapes) +
    scale_alpha_manual(values = point.alpha, guide="none") +
    scale_size_manual(values = point.sizes, guide="none") +
    scale_color_manual(values = r2_colors) +
    scale_fill_manual(values = r2_colors, guide="none") +
    xlim(loc$xrange/1e6) +
    labs(x = "", y = "PIP", shape = "Type", color = expression(R^2)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, l=10, t=0, r=10))

  if (!is.null(highlight_pip)) {
    p_pip <- p_pip +
      geom_hline(yintercept=highlight_pip, linetype="dashed", color = "red")
  }

  # QTL panel
  p_qtl <- ggplot(finemap_qtl_res, aes(x=.data$pos/1e6)) +
    geom_rect(aes(xmin=loc$xrange[1]/1e6, xmax=loc$xrange[2]/1e6, ymin=0, ymax=1),
              fill="gray90") +
    geom_segment(aes(x=.data$pos/1e6, xend=.data$pos/1e6, y=0, yend=1), color="salmon") +
    labs(title = paste(focus_gene_name, focus_gene_context, focus_gene_type),
         y = "QTL") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size=10),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      # strip.text.y.left = element_text(angle=0, size=8),
      # strip.background = element_blank(),
      plot.margin = margin(b=0, l=10, t=0, r=10))

  if (label_QTLs){
    p_qtl <- p_qtl +
      geom_text_repel(aes(y=0.5, label=.data$id), size=label.text.size, color="black")
  }

  # gene track panel
  p_genes <- gg_genetracks(loc, xticks=FALSE,
                           filter_gene_biotype = filter_genetrack_biotype,
                           text_pos="top") +
    labs(x = paste0("chr", chrom)) +
    theme_bw() +
    theme(legend.position = "none",
          panel.border= element_blank(),
          plot.margin = margin(b=10, l=10, t=0, r=10))

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
      geom_vline(xintercept = highlight_pos/1e6, linetype="dotted", color = "blue", size=0.5)
  }

  cowplot::plot_grid(p_pvalue, p_pip, p_qtl, p_genes, ncol = 1,
                     rel_heights = panel.heights, align = "v", axis="tblr")
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

  # estimated group PVE (all iterations)
  group_pve_iters <- group_prior_var_iters*group_prior_iters*group_size/gwas_n
  # estimated enrichment of genes (all iterations)
  enrichment_iters <- t(sapply(rownames(group_prior_iters)[rownames(group_prior_iters)!="SNP"], function(x){
    group_prior_iters[rownames(group_prior_iters)==x,]/group_prior_iters[rownames(group_prior_iters)=="SNP"]}))

  # make convergence plots

  # inclusion plot
  df <- data.frame(niter = rep(1:ncol(group_prior_iters), nrow(group_prior_iters)),
                   value = unlist(lapply(1:nrow(group_prior_iters), function(x){group_prior_iters[x,]})),
                   group = rep(rownames(group_prior_iters), each=ncol(group_prior_iters)))
  factor_levels <- c(setdiff(rownames(group_prior_iters), "SNP"), "SNP")
  df$group <- factor(df$group, levels = factor_levels)

  p_pi <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab(bquote(pi)) +
    ggtitle("Proportion Causal") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "Group")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size))

  # effect size plot
  df <- data.frame(niter = rep(1:ncol(group_prior_var_iters), nrow(group_prior_var_iters)),
                   value = unlist(lapply(1:nrow(group_prior_var_iters), function(x){group_prior_var_iters[x,]})),
                   group = rep(rownames(group_prior_var_iters), each=ncol(group_prior_var_iters)))
  df$group <- factor(df$group, levels = factor_levels)

  p_sigma2 <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab(bquote(sigma^2)) +
    ggtitle("Effect Size") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "Group")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size))

  # PVE plot
  df <- data.frame(niter = rep(1:ncol(group_pve_iters), nrow(group_pve_iters)),
                   value = unlist(lapply(1:nrow(group_pve_iters), function(x){group_pve_iters[x,]})),
                   group = rep(rownames(group_pve_iters), each=ncol(group_pve_iters)))
  df$group <- factor(df$group, levels = factor_levels)

  p_pve <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab(bquote(h[G]^2)) +
    ggtitle("PVE") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "Group")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size))

  # enrichment plot
  df <- data.frame(niter = rep(1:ncol(enrichment_iters), nrow(enrichment_iters)),
                   value = unlist(lapply(1:nrow(enrichment_iters), function(x){enrichment_iters[x,]})),
                   group = rep(rownames(enrichment_iters), each=ncol(enrichment_iters)))
  df$group <- factor(df$group, levels = factor_levels[factor_levels!="SNP"])

  p_enrich <- ggplot(df, aes(x=.data$niter, y=.data$value, group=.data$group, color=.data$group)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colors) +
    xlab("Iteration") + ylab(bquote(pi[G]/pi[S])) +
    ggtitle("Enrichment") +
    theme_cowplot() +
    theme(plot.title=element_text(size=title.size)) +
    expand_limits(y=0) +
    guides(color = guide_legend(title = "Group")) +
    theme(legend.title = element_text(size=legend.size, face="bold"),
          legend.text = element_text(size=legend.size))

  cowplot::plot_grid(p_pi, p_sigma2, p_enrich, p_pve, ncol = ncol)
}
