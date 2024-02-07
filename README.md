# cTWAS: integrating molecular QTLs and GWAS for gene discovery

Expression Quantitative Trait Loci (eQTLs) have often been used to nominate candidate genes from Genome-wide association studies (GWAS). However, commonly used methods are susceptible to false positives largely due to Linkage Disequilibrium of eQTLs with causal variants acting on the phenotype directly. Our method, causal-TWAS (cTWAS), addressed this challenge by borrowing ideas from statistical fine-mapping. It is a generalization of Transcriptome-wide association studies (TWAS), but when analyzing any gene, it adjusts for other nearby genes and all nearby genetic variants.  

Running cTWAS involves four main steps: preparing input data, imputing gene z-scores, estimating parameters, and fine-mapping genes and variants. The output of cTWAS are posterior inclusion probabilities (PIPs) for all variants and genes with expression models. We have included a [tutorial](https://xinhe-lab.github.io/ctwas/articles/ctwas_summary_statistics.html) of how to use the ctwas software. You can [browse source code](https://github.com/xinhe-lab/ctwas) and [report a bug](https://github.com/xinhe-lab/ctwas/issues) here. 

<img style="display:block;margin:auto" width="700" height="300" src="./workflow.png">

## Install

Install `ctwas`:

```
remotes::install_github("xinhe-lab/ctwas",ref = "main")
```

Currently, `ctwas` has only been tested on linux systems. 
We recommend running `ctwas` on a High-Performance Computing system.

## Citing this work

If you find the `ctwas` package or any of the source code in this
repository useful for your work, please cite:

> Zhao S, Crouse W, Qian S, Luo K, Stephens M, He X. 
> Adjusting for genetic confounders in transcriptome-wide association 
> studies improves discovery of risk genes of complex traits. 
> Nat Genet (2024). https://doi.org/10.1038/s41588-023-01648-9

## Useful resources

We have pre-computed the LD matrices of European samples from UK Biobank. 
They can be downloaded [here](https://uchicago.box.com/s/jqocacd2fulskmhoqnasrknbt59x3xkn). 

cTWAS requires the expression prediction models, or weights, of genes. 
The pre-computed weights of GTEx expression and splicing traits can be downloaded from [PredictDB](https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/). 

## Acknowledgement

We acknowledge the authors of `susieR` package for using their codes.

Original `susieR` code obtained by:
```
git clone git@github.com:stephenslab/susieR.git
git checkout c7934c0
```

Minor edits to make it accept different prior variances for each variable.
