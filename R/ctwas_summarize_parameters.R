#' Summarize and plot cTWAS parameter estimates
#'
#' @param param a list of cTWAS parameter estimation result from \code{est_param}
#'
#' @param gwas_n the sample size of the GWAS summary statistics
#'
#' @param plot if TRUE, return a plot of the estimated parameters
#'
#' @param title.size font size of the plot title
#'
#' @param legend.size font size of the plot legend title
#'
#' @importFrom logging loginfo
#' @import ggplot2
#' @importFrom cowplot plot_grid theme_cowplot
#'
#' @export
#'
summarize_param <- function(param,
                            gwas_n = NA,
                            plot = TRUE,
                            title.size = 10,
                            legend.size = 8){

  group_prior <- param$group_prior
  # estimated group prior (all iterations)
  group_prior_iters <- param$group_prior_iters

  group_prior_var <- param$group_prior_var
  # estimated group prior variance (all iterations)
  group_prior_var_iters <- param$group_prior_var_iters

  # set group size
  group_size <- param$group_size
  group_size <- group_size[rownames(group_prior_iters)]

  # estimated group PVE (all iterations)
  group_pve_iters <- group_prior_var_iters*group_prior_iters*group_size/gwas_n
  group_pve <- group_pve_iters[,ncol(group_pve_iters)]

  # estimated enrichment of genes (all iterations)
  enrichment_iters <- t(sapply(rownames(group_prior_iters)[rownames(group_prior_iters)!="SNP"], function(x){
    group_prior_iters[rownames(group_prior_iters)==x,]/group_prior_iters[rownames(group_prior_iters)=="SNP"]}))
  enrichment <- enrichment_iters[,ncol(enrichment_iters)]

  outlist <- list(group_size = group_size,
                  group_prior = group_prior,
                  group_prior_var = group_prior_var,
                  enrichment = enrichment,
                  group_pve = group_pve,
                  total_pve = sum(group_pve),
                  attributable_pve = group_pve/sum(group_pve))

  if (plot){

    # inclusion plot
    df <- data.frame(niter = rep(1:ncol(group_prior_iters), nrow(group_prior_iters)),
                     value = unlist(lapply(1:nrow(group_prior_iters), function(x){group_prior_iters[x,]})),
                     group = rep(rownames(group_prior_iters), each=ncol(group_prior_iters)))
    factor_levels <- c(setdiff(rownames(group_prior_iters), "SNP"), "SNP")
    df$group <- factor(df$group, levels = factor_levels)

    p_pi <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
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

    p_sigma2 <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
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

    p_pve <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
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

    p_enrich <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(pi[G]/pi[S])) +
      ggtitle("Enrichment") +
      theme_cowplot() +
      theme(plot.title=element_text(size=title.size)) +
      expand_limits(y=0) +
      guides(color = guide_legend(title = "Group")) +
      theme(legend.title = element_text(size=legend.size, face="bold"),
            legend.text = element_text(size=legend.size)) +
      scale_colour_discrete(drop = FALSE)

    outlist$convergence_plot <- cowplot::plot_grid(p_pi, p_sigma2, p_enrich, p_pve)
  }

  return(outlist)
}
