# cTWAS

The ctwas package implements a novel statistical framework for integrating eQTL and GWAS data, enabling reliable gene discoveries. Our approach can be viewed as a generalization of TWAS, which we term "causal-TWAS" (cTWAS).

Running cTWAS involves four main steps: preparing input data, imputing gene z-scores, estimating parameters, and fine-mapping the genes and variants. The output of cTWAS is a posterior inclusion probability (PIP) for each variant and each gene with an expression model. We have included a [tutorial](https://xinhe-lab.github.io/ctwas/articles/transition.html) of how to use the ctwas software. 

For more details, Please read the [manuscript](https://doi.org/10.1101/2022.09.27.509700) on bioRxiv.

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

## Acknowledgement

We acknowledge the authors of `susieR` package for using their codes.

Original `susieR` code obtained by:
```
git clone git@github.com:stephenslab/susieR.git
git checkout c7934c0
```

Minor edits to make it accept different prior variances for each variable.