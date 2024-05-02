
#' Plot the cTWAS result for a single locus
#'
#' @param finemap_res a data frame of cTWAS finemapping result
#' @param region_id region ID to be plotted
#' @param weights a list of weights
#' @param screened_region_data a list of screened region_data for regions with strong signals
#' @param cor_matrix_path path of correlation matrices
#' @param ens_db Ensembl database
#' @param locus_range a vector of start and end positions to define the locus region to be plotted
#' If NULL, the entire region will be plotted (with 100 bp flanks on both sides).
#' @param focus_gene the gene to focus the plot. If NULL, use the lead TWAS feature
#' @param gene_biotype specified biotypes to be displayed in gene tracks or NULL to display all possible ones.
#' By default, only show protein coding genes.
#' @param use_cowplot TRUE/FALSE. If TRUE, use cowplot::plot_grid to align the panels.
#' Otherwise, use gtable to align the panels.
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import gtable
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
                           use_cowplot = FALSE) {

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

  # add the correlations to the focus gene
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

  # fix the colors: 0.4~1: purple, lead gene: salmon, others "#7FC97F"
  # may be modified to add more colors representing different r2
  finemap_region_res$r2_category <- cut(finemap_region_res$r2,
                                        breaks = c(0, 0.4, 1, Inf),
                                        labels = c("0-0.4", "0.4-1", "1"))
  finemap_region_res$r2_color <- factor(finemap_region_res$r2_category,
                                        levels = c("0-0.4", "0.4-1", "1"),
                                        labels = c("#7FC97F", "purple", "salmon"))

  # create a locus object for plotting
  loc <- locuszoomr::locus(
    data = finemap_region_res,
    seqname = unique(finemap_region_res$chrom),
    xrange = locus_range,
    ens_db = ens_db,
    labs = "id")

  loc$data$r2_color <- finemap_region_res$r2_color

  # p-value panel
  p_pvalue <- ggplot(loc$data, aes(x=pos/1e6, y=p,  #### modified: loc$data
                                   shape=type,
                                   label=gene_name,
                                   size=object_type,
                                   alpha=object_type,
                                   fill=r2_color)) +   # Use the new color mapping variable
    geom_point() +
    geom_point(aes(x=pos/1e6, y=-log10(p)), color="black") +
    geom_text_repel(size=2) +
    scale_fill_identity() +  # Use the colors directly
    scale_alpha_manual(values=c("SNP"=0.7, "non-SNP"=1), guide="none") +
    scale_size_manual(values=c("SNP"=1, "non-SNP"=3.5), guide="none") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    xlim(loc$xrange/1e6) +
    xlab("") +
    ylab("-log10(p-value)") +
    theme_bw() +
    theme(legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, t=0))

  # PIP panel
  p_pip <- ggplot(loc$data, aes(x=pos/1e6, y=susie_pip,  #### modified: loc$data
                                shape=type,
                                label=gene_name,
                                size=object_type,
                                alpha=object_type,
                                fill=r2_color)) +  # Use the new color mapping variable
    geom_point() +
    geom_point(aes(x=pos/1e6, y=susie_pip), color="black") +
    geom_text_repel(size=2) +
    scale_fill_identity() +  # Use the colors directly from the r2_color variable
    scale_alpha_manual(values=c("SNP"=0.7, "non-SNP"=1), guide="none") +
    scale_size_manual(values=c("SNP"=1, "non-SNP"=3.5), guide="none") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    xlim(loc$xrange/1e6) +
    xlab("") +
    ylab("SuSiE PIP") +
    theme_bw() +
    theme(legend.position = "none",
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, t=0))

  # QTL panel
  # get QTL info
  snps <- rownames(weights[[focus_gene]]$wgt)
  snps <- finemap_region_res[finemap_region_res$id %in% snps, ]
  snps$qtl_name <- NULL

  p_qtl <- ggplot(snps) +
    geom_rect(aes(xmin=locus_range[1]/1e6, xmax=locus_range[2]/1e6, ymin=0.1, ymax=0.2),
              fill = "lightgray") +
    geom_segment(aes(x=pos/1e6, xend=pos/1e6, y=0.1, yend=0.2)) +
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
      plot.margin = margin(b=0, l=-6, t=0))

  # gene track panel
  p_genes <- gg_genetracks(loc, xticks=FALSE, filter_gene_biotype = gene_biotype, text_pos="top") #+ theme(plot.margin=margin(t=-5))

  if (use_cowplot == TRUE){
    cowplot::plot_grid(p_pvalue, p_pip, p_qtl, p_genes, ncol = 1,
                       rel_heights = c(4, 4, 0.3, 4), align = "v")
  } else {
    # align panels
    # https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
    g_pvalue <- ggplotGrob(p_pvalue)
    g_pip <- ggplotGrob(p_pip)
    g_qtl <- ggplotGrob(p_qtl)
    g_genes <- ggplotGrob(p_genes)

    # bind the gtables
    g <- rbind(g_pvalue, g_pip, size="max")
    panels_extent <- g %>% find_panel()

    pg <- g %>%
      gtable_add_rows(heights=unit(0.2, "cm"), pos = 0) %>% ## modified: Add space at the very top
      gtable_add_rows(heights=unit(0.5, "cm")) %>%
      gtable_add_grob(g_qtl, t=-1, b=-1, r=panels_extent$r-2, l=panels_extent$l) %>%
      gtable_add_rows(heights = unit(0.1, "cm")) %>%
      gtable_add_rows(heights = unit(5, "cm")) %>%
      gtable_add_grob(g_genes, t=-1, b=-1, r=panels_extent$r+1, l=panels_extent$l)

    pg$widths <- unit.pmax(g_pvalue$widths, g_pip$widths)
    grid.newpage()
    grid.draw(pg)
  }
}
