# cTWAS: integrating molecular QTLs and GWAS for gene discovery

Expression Quantitative Trait Loci (eQTLs) have often been used to nominate candidate genes from Genome-wide association studies (GWAS). However, commonly used methods are susceptiable to false positives largely due to Linkage Disequilibrium of eQTLs with causal variants acting on the phenotype directly. Our method, causal-TWAS (cTWAS), addressed this challenge by borrowing ideas from statistical fine-mapping. It is a generalization of Transcriptome-wide association studies (TWAS), but when analyzing any gene, it adjusts for other nearby genes and all nearby geentic variants.  

Running cTWAS involves four main steps: preparing input data, imputing gene z-scores, estimating parameters, and fine-mapping genes and variants. The output of cTWAS are posterior inclusion probabilities (PIPs) for all variants and genes with expression models. We have included a [tutorial](https://xinhe-lab.github.io/ctwas/articles/transition.html) of how to use the ctwas software. 

![](./workflow.png =250x)

## Install

Install `ctwas`:

```
remotes::install_github("xinhe-lab/ctwas",ref = "main")
```
Currently, ctwas can only be installed on linux systems. We suggset to run ctwas using a high throughput computing system.

## Citing this work

If you find the `ctwas` package or any of the source code in this
repository useful for your work, please cite:

> Siming, Z., Wesley, C., Sheng, Q., Kaixuan, L., Stephens, M. & Xin, H. (2022). 
> Adjusting for genetic confounders in transcriptome-wide association 
> studies leads to reliable detection of causal genes
> application to genetic fine mapping. *bioRxiv
> https://doi.org/10.1101/2022.09.27.509700

## Useful resources

We have pre-computed the LD matrices of European samples from UK Biobank. They can be downloaded [here](https://uchicago.box.com/s/jqocacd2fulskmhoqnasrknbt59x3xkn). 

cTWAS requires the expression prediction models, or weights, of genes. The pre-computed weights of GTEx expression and splicing traits can be downloaded from [PredictDB](https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/). 

To run cTWAS, it is useful to perform some pre-processing on the GWAS summary statistics, to check allele flipping and filter problemetic variants due to LD mismatch between reference and in-sample LD. Some useful software for this purpose, include [SuSiE-RSS](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html) and [DENTIST](https://github.com/Yves-CHEN/DENTIST/). 

## Acknowledgement

We acknowledge the authors of `susieR` package for using their codes.

Original `susieR` code obtained by:
```
git clone git@github.com:stephenslab/susieR.git
git checkout c7934c0
```

Minor edits to make it accept different prior variances for each variable.