#' Summarize and plot cTWAS parameter estimates
#' 
#' @param outputdir a string, the directory where ctwas output is stored
#' 
#' @param outname a string, the output name
#' 
#' @param gwas_n the sample size of the GWAS summary statistics
#' 
#' @param thin the proportion of SNPs used for parameter estimation, as supplied to \code{ctwas_rss}
#' 
#' @param plot_estimates if TRUE, return a plot of the estimated parameters
#'
#' @export
#' 
ctwas_summarize_parameters <- function(outputdir,
                                       outname,
                                       gwas_n = NA,
                                       thin = 1,
                                       plot_estimates = T){
  library(ggplot2)
  library(cowplot)
  
  load(paste0(outputdir, outname, ".s2.susieIrssres.Rd"))
  
  #estimated group prior (all iterations)
  estimated_group_prior_all <- group_prior_rec
  estimated_group_prior_all["SNP",] <- estimated_group_prior_all["SNP",]*thin #adjust parameter to account for thin argument
  group_prior <- estimated_group_prior_all[,ncol(estimated_group_prior_all)]
  
  #estimated group prior variance (all iterations)
  estimated_group_prior_var_all <- group_prior_var_rec
  group_prior_var <- estimated_group_prior_var_all[,ncol(estimated_group_prior_var_all)]
  
  #set group size
  ctwas_res_s1 <- as.data.frame(data.table::fread(paste0(outputdir, outname, ".s1.susieIrss.txt"), header = T))
  
  group_size <- table(ctwas_res_s1$type)
  group_size["SNP"] <- group_size["SNP"]/thin #adjust to account for thin argument
  
  group_size <- group_size[rownames(estimated_group_prior_all)]
  
  #estimated group PVE (all iterations)
  estimated_group_pve_all <- estimated_group_prior_var_all*estimated_group_prior_all*group_size/gwas_n
  group_pve <- estimated_group_pve_all[,ncol(estimated_group_pve_all)]
  
  #estimated enrichment of genes (all iterations)
  estimated_enrichment_all <- t(sapply(rownames(estimated_group_prior_all)[rownames(estimated_group_prior_all)!="SNP"], 
                                       function(x){estimated_group_prior_all[rownames(estimated_group_prior_all)==x,]/estimated_group_prior_all[rownames(estimated_group_prior_all)=="SNP"]}))
  enrichment <- estimated_enrichment_all[,ncol(estimated_enrichment_all)]
  
  outlist <- list(group_size = group_size,
                  group_prior = group_prior,
                  group_prior_var = group_prior_var,
                  enrichment = enrichment,
                  group_pve = group_pve,
                  total_pve = sum(group_pve),
                  attributable_pve = group_pve/sum(group_pve))
  
  if (plot_estimates){
    title_size <- 12
    
    #inclusion plot
    df <- data.frame(niter = rep(1:ncol(estimated_group_prior_all), nrow(estimated_group_prior_all)),
                     value = unlist(lapply(1:nrow(estimated_group_prior_all), function(x){estimated_group_prior_all[x,]})),
                     group = rep(rownames(estimated_group_prior_all), each=ncol(estimated_group_prior_all)))
    df$group <- as.factor(df$group)
    factor_levels <- levels(df$group)
    
    p_pi <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(pi)) +
      ggtitle("Proportion Causal") +
      theme_cowplot()
    
    p_pi <- p_pi + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold"))
    
    #effect size plot
    df <- data.frame(niter = rep(1:ncol(estimated_group_prior_var_all), nrow(estimated_group_prior_var_all)),
                     value = unlist(lapply(1:nrow(estimated_group_prior_var_all), function(x){estimated_group_prior_var_all[x,]})),
                     group = rep(rownames(estimated_group_prior_var_all), each=ncol(estimated_group_prior_var_all)))
    df$group <- as.factor(df$group)
    
    p_sigma2 <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(sigma^2)) +
      ggtitle("Effect Size") +
      theme_cowplot()
    
    p_sigma2 <- p_sigma2 + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold"))
    
    #PVE plot
    df <- data.frame(niter = rep(1:ncol(estimated_group_pve_all), nrow(estimated_group_pve_all)),
                     value = unlist(lapply(1:nrow(estimated_group_pve_all), function(x){estimated_group_pve_all[x,]})),
                     group = rep(rownames(estimated_group_pve_all), each=ncol(estimated_group_pve_all)))
    df$group <- as.factor(df$group)
    
    p_pve <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(h[G]^2)) +
      ggtitle("PVE") +
      theme_cowplot()
    
    p_pve <- p_pve + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold"))
    
    #enrichment plot
    df <- data.frame(niter = rep(1:ncol(estimated_enrichment_all), nrow(estimated_enrichment_all)),
                     value = unlist(lapply(1:nrow(estimated_enrichment_all), function(x){estimated_enrichment_all[x,]})),
                     group = rep(rownames(estimated_enrichment_all), each=ncol(estimated_enrichment_all)))
    df$group <- as.factor(df$group)
    #levels(df$group) <- factor_levels
    
    p_enrich <- ggplot(df, aes(x=niter, y=value, group=group)) +
      geom_line(aes(color=group)) +
      geom_point(aes(color=group)) +
      xlab("Iteration") + ylab(bquote(pi[G]/pi[S])) +
      ggtitle("Enrichment") +
      theme_cowplot()
    
    p_enrich <- p_enrich + theme(plot.title=element_text(size=title_size)) + 
      expand_limits(y=0) + 
      guides(color = guide_legend(title = "Group")) + theme (legend.title = element_text(size=12, face="bold")) + scale_colour_discrete(drop = FALSE)
    
    convergence_plot <- plot_grid(p_pi, p_sigma2, p_enrich, p_pve)
    
    outlist[["convergence_plot"]] <- convergence_plot
  }
  
  return(outlist)
}
