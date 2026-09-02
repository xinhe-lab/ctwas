# Plots the cTWAS result for a single locus

Plots the cTWAS result for a single locus

## Usage

``` r
make_locusplot(
  finemap_res,
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
  LD.breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  LD.colors = c("grey", "cyan", "green", "orange", "red", "salmon"),
  L = 5,
  cs.colors = c("firebrick", "dodgerblue", "forestgreen", "darkmagenta", "darkorange",
    "grey"),
  focal.colors = c("grey", "salmon"),
  label_QTLs = TRUE,
  highlight_pval = NULL,
  highlight_pip = 0.8,
  highlight_pos = NULL,
  highlight.color = "red",
  highlight.pos.color = "magenta",
  point.sizes = c(1, 3),
  point.alpha = c(0.4, 0.6),
  point.shapes = c(16, 15, 18, 17, 10, 12, 14, 11),
  label.text.size = 2.5,
  max.overlaps = 10,
  axis.text.size = 8,
  axis.title.size = 10,
  legend.text.size = 9,
  legend.position = "top",
  legend.nrow = c(1, 1),
  genelabel.cex.text = 0.7,
  panel.heights = c(4, 4, 1, 4),
  verbose = FALSE
)
```

## Arguments

  - finemap\_res:
    
    a data frame of cTWAS finemapping result

  - region\_id:
    
    region ID to be plotted

  - ens\_db:
    
    Ensembl database

  - weights:
    
    a list of proprocessed weights, used to plot QTL track. If NULL,
    will not plot the QTL track.

  - R\_snp\_gene:
    
    SNP-gene correlation matrix of the region. If both R\_snp\_gene and
    R\_gene are available, color data points with correlations with the
    focal gene.

  - R\_gene:
    
    gene-gene correlation matrix of the region

  - locus\_range:
    
    a vector of start and end positions to define the region boundary to
    be plotted. If NULL, the entire region will be plotted (with 100 bp
    flanks on both sides).

  - focal\_id:
    
    the focal ID.

  - focal\_gene:
    
    the focal gene name. By default, use the gene with the highest PIP.

  - filter\_protein\_coding\_genes:
    
    If TRUE, limits to protein coding genes only.

  - filter\_cs:
    
    If TRUE, limits to credible sets.

  - color\_pval\_by:
    
    Options to color the p-value track. "LD": colors the p-value track
    by the correlations to the focal gene. "cs": colors the p-value
    track by the credible sets. "none": uses the same color for
    non-focal genes.

  - color\_pip\_by:
    
    Options to color the PIP track. "LD": colors the PIP track by the
    correlations to the focal gene. "cs": colors the PIP track by the
    credible sets. "none": uses the same color for non-focal genes.

  - LD.breaks:
    
    Breaks of LD intervals.

  - LD.colors:
    
    Colors for correlation levels.

  - L:
    
    Number of effects in finemapping.

  - cs.colors:
    
    Colors for credible sets.

  - focal.colors:
    
    Colors for non-focal and focal gene.

  - label\_QTLs:
    
    If TRUE, label SNP IDs in the QTL panel.

  - highlight\_pval:
    
    p-value to highlight with a horizontal line

  - highlight\_pip:
    
    PIP to highlight with a horizontal line

  - highlight\_pos:
    
    genomic positions to highlight with vertical lines

  - highlight.color:
    
    color for horizontal lines highlighting p-values and PIPs.

  - highlight.pos.color:
    
    color for vertical lines highlighting genomic positions.

  - point.sizes:
    
    size values for SNP and non-SNP data points in the scatter plots

  - point.alpha:
    
    alpha values for SNP and non-SNP data points in the scatter plots

  - point.shapes:
    
    shapes values for data points of different types in the scatter
    plots

  - label.text.size:
    
    Font size for gene and SNP label text

  - max.overlaps:
    
    Setting for geom\_text\_repel() function to label texts. Exclude
    text labels when they overlap too many other things.

  - axis.text.size:
    
    Font size for axis label text.

  - axis.title.size:
    
    Font size for axis title text.

  - legend.text.size:
    
    Font size for legend text.

  - legend.position:
    
    position to put legends. If "none", no legends will be shown.

  - legend.nrow:
    
    Number of rows for legend text in the panels.

  - genelabel.cex.text:
    
    Controls the size of the gene labels.

  - panel.heights:
    
    Relative heights of the panels.

  - verbose:
    
    If TRUE, print detail messages.
