# cTWAS: integrating molecular QTLs and GWAS for gene discovery

Expression Quantitative Trait Loci (eQTLs) have often been used to nominate candidate genes from Genome-wide association studies (GWAS). However, commonly used methods are susceptible to false positives largely due to Linkage Disequilibrium of eQTLs with causal variants acting on the phenotype directly. Our method, causal-TWAS (cTWAS), addressed this challenge by borrowing ideas from statistical fine-mapping. It is a generalization of Transcriptome-wide association studies (TWAS), but when analyzing any gene, it adjusts for other nearby genes and all nearby genetic variants.  

## New version: multi-group cTWAS

We have updated the `ctwas` package. The new version, multi-group cTWAS (or multi-cTWAS), can analyze GWAS and molecular QTL data from multiple modalities and contexts. 
You can follow the new [tutorials](https://xinhe-lab.github.io/multigroup_ctwas/). (Tutorial for the old version can be found [here](https://xinhe-lab.github.io/ctwas/articles/ctwas_summary_statistics.html)).

You can browse [source code for the new version](https://github.com/xinhe-lab/ctwas/tree/multigroup) and report a bug [here](https://github.com/xinhe-lab/ctwas/issues). You can also join our [Google Group](https://groups.google.com/g/ctwas_users) to receive notifications of software updates. 

## Citing this work

If you find the `ctwas` package or any of the source code in this
repository useful for your work, please cite:

> Zhao S, Crouse W, Qian S, Luo K, Stephens M, He X. 
> Adjusting for genetic confounders in transcriptome-wide association 
> studies improves discovery of risk genes of complex traits. 
> Nat Genet (2024). https://doi.org/10.1038/s41588-023-01648-9

