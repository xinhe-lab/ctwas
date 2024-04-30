
#' Plot the output of cTWAS at a single locus
#'
#' @export
#' 
ctwas_locusplot <- function(finemap_res,
                            screened_region_data,
                            region_id,
                            cor_matrix_path = getwd(),
                            locus_range = NULL,
                            focus_gene = NULL,
                            gene_biotype="protein_coding",
                            genome_build = "hg38") {

  finemap_res <- finemap_res[finemap_res$region_id==region_id,] #select target region
  finemap_res$p <- (-log(2) - pnorm(abs(finemap_res$z), lower.tail=F, log.p=T))/log(10) #add pvalue
  finemap_res$object_type <- finemap_res$type #adding object_type for size and alpha. SNP and non-SNP
  finemap_res$object_type[finemap_res$object_type!="SNP"] <- "non-SNP"

  if (!is.null(focus_gene)) { #if not specified, selecting the top gene
    focus_gene <- finemap_res[which(finemap_res$gene_name == focus_gene), ]$id
  } else{
    focus_gene <- finemap_res$id[which(finemap_res$susie_pip == max(finemap_res$susie_pip[finemap_res$type != "SNP"]))[1]]
  }

  if (!is.null(locus_range)){ #if not specified, plotting the whole region
    finemap_res <- finemap_res[finemap_res$pos>=locus_range[1] & finemap_res$pos<=locus_range[2],, drop=F]
  }else{
    locus_range <- c(min(finemap_res$pos)-1000, max(finemap_res$pos))
  }

  ###Adding the correlation info, with respect to the top gene
  R_gene <- readRDS(paste0(cor_matrix_path,"region.",region_id,".R_gene.RDS"))
  R_snp_gene <- readRDS(paste0(cor_matrix_path,"region.",region_id,".R_snp_gene.RDS"))
  R_snp <-  readRDS(paste0(cor_matrix_path,"region.",region_id,".R_snp.RDS"))

  rownames(R_gene) <- screened_region_data[[region_id]]$gid
  colnames(R_gene) <- screened_region_data[[region_id]]$gid
  rownames(R_snp_gene) <- screened_region_data[[region_id]]$sid
  colnames(R_snp_gene) <-screened_region_data[[region_id]]$gid
  rownames(R_snp) <- screened_region_data[[region_id]]$sid
  colnames(R_snp) <- screened_region_data[[region_id]]$sid

  finemap_res$r2max <- NA
  finemap_res$r2max[finemap_res$type!="SNP"] <- R_gene[focus_gene, finemap_res$id[finemap_res$type!="SNP"]]
  finemap_res$r2max[finemap_res$type=="SNP"] <- R_snp_gene[finemap_res$id[finemap_res$type=="SNP"], focus_gene]

  r2cut <- 0.4
  finemap_res$r2max<-finemap_res$r2max > r2cut
  finemap_res$r2max[finemap_res$type != "SNP" & finemap_res$id == focus_gene] <- 2

  #loading gene info database
  genome_build = "hg38"
  if(genome_build=="hg38"){
    library(EnsDb.Hsapiens.v86)
    ensdb <- EnsDb.Hsapiens.v86
  }else if(genome_build=="hg19"){
    library(EnsDb.Hsapiens.v75)
    ensdb <- EnsDb.Hsapiens.v75
  }

  # create a locus object for plotting
  loc <- locus(
    data = finemap_res,
    seqname = unique(finemap_res$chrom),
    xrange = locus_range,
    ens_db = ensdb,
    labs = "id"
  )

  ### modified:
  # fix the colors: 0.4~1: purple, lead gene: salmon, others "#7FC97F"
  # may be modified to add more colors representing different ld r2
  loc$data$r2max_color <- cut(loc$data$r2max,
                              breaks = c(-Inf, 0.4, 1, Inf),
                              labels = c("#7FC97F", "purple", "salmon"))


  p1 <- ggplot(loc$data, aes(x=pos/1e6, y=p,  #### modified: loc$data
                             shape=type,
                             label=gene_name,
                             size=object_type,
                             alpha=object_type,
                             fill=r2max_color)) +   # Use the new color mapping variable
    geom_point() +
    geom_point(aes(x=pos/1e6, y=-log10(p)), color="black") +
    geom_text_repel(size=2) +
    scale_fill_identity() +  # Use the colors directly
    scale_alpha_manual(values=c("SNP"=0.7, "non-SNP"=1), guide="none") +
    scale_size_manual(values=c("SNP"=1, "non-SNP"=3.5), guide="none") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    ylab("-log10(p-value)") + xlim(loc$xrange/1e6) +
    xlab("") + theme_bw() + theme(legend.position = "none",
                                  axis.ticks.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  panel.border= element_blank(),
                                  axis.line = element_line(colour = "black"),
                                  plot.margin = margin(b=0, t=0))

  p2 <- ggplot(loc$data, aes(x=pos/1e6, y=susie_pip,  #### modified: loc$data
                             shape=type,
                             label=gene_name,
                             size=object_type,
                             alpha=object_type,
                             fill=r2max_color)) +  # Use the new color mapping variable
    geom_point() +
    geom_point(aes(x=pos/1e6, y=susie_pip), color="black") +
    geom_text_repel(size=2) +
    scale_fill_identity() +  # Use the colors directly from the r2max_color variable
    scale_alpha_manual(values=c("SNP"=0.7, "non-SNP"=1), guide="none") +
    scale_size_manual(values=c("SNP"=1, "non-SNP"=3.5), guide="none") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    ylab("SuSiE PIP") + xlim(loc$xrange/1e6) +
    xlab("") + theme_bw() +
    theme(legend.position = "none",
          panel.border= element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(b=0, t=0))

  ## get qtl info
  snps <- rownames(weights[[focus_gene]]$wgt)
  snps <- finemap_res[finemap_res$id %in% snps, ]
  snps$qtl_name <- NULL

  ## draw qtl tracks
  x_start=locus_range[1]/1e6
  x_end=locus_range[2]/1e6

  q1<- ggplot(snps) +
    geom_rect(mapping=aes(xmin=x_start, xmax=x_end,
                          ymin=0.1, ymax=0.2)) +
    geom_segment(mapping=aes(x=pos/1e6,
                             xend=pos/1e6,
                             y=0.1, yend=0.2)) +
    theme(
      axis.title = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.text.y.left = element_text(
        angle=0,
        size=8),
      strip.background = element_blank(),
      plot.margin = margin(
        b=0,
        l=-6,
        t=0)
    )

  g4 <- gg_genetracks(loc, xticks=FALSE, filter_gene_biotype = gene_biotype, text_pos="top") #+ theme(plot.margin=margin(t=-5))

  ## Adding gene tracks
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g3 <- ggplotGrob(q1)
  g<-rbind(g1, g2, size="max")
  panels_extent <- g %>% find_panel()

  pg <- g %>%
    gtable_add_rows(heights=unit(0.2, "cm"), pos = 0) %>% # ## modified: Add space at the very top
    gtable_add_rows(heights=unit(0.5, "cm")) %>%
    gtable_add_grob(g3,
                    t=-1, b=-1,
                    r=panels_extent$r-2, l=panels_extent$l) %>%
    gtable_add_rows(heights = unit(0.1, "cm")) %>%
    gtable_add_rows(heights = unit(5, "cm")) %>%
    gtable_add_grob(ggplotGrob(g4),
                    t=-1, b=-1,
                    r=panels_extent$r+1, l=panels_extent$l)  ## modified: I always have null g4

  pg$widths <- unit.pmax(g2$widths, g1$widths)
  grid.newpage() 
  grid.draw(pg)
}