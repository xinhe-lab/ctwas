% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctwas_convert_geno_to_LD_matrix.R
\name{convert_geno_to_LD_matrix}
\alias{convert_geno_to_LD_matrix}
\title{Converts PLINK genotype data to LD matrices and SNP info files,
saves LD matrices as .RDS files and SNP info as .Rvar files}
\usage{
convert_geno_to_LD_matrix(
  region_info,
  genotype_files,
  varinfo_files,
  chrom = 1:22,
  outputdir = getwd(),
  outname = "",
  include_variance = TRUE,
  include_allele_freq = TRUE,
  show_progress_bar = TRUE,
  verbose = FALSE,
  logfile = NULL
)
}
\arguments{
\item{region_info}{a data frame of region definitions, with columns:
chrom, start, stop, and region_id.}

\item{genotype_files}{Reference genotype files in PLINK binary genotype data
in .pgen or .bed format. It should contain files for all chromosomes
(from 1 to 22), one file per chromosome.}

\item{varinfo_files}{Reference variant information files in PLINK
.pvar or .bim format. It could have one file per chromosome or
have one big file for all chromosomes.
The output will use the genome positions in \code{varinfo_files}.}

\item{chrom}{a vector of chromosome numbers to process genotype data.}

\item{outputdir}{Output directory.}

\item{outname}{Output filestem.}

\item{include_variance}{If TRUE, include variance in .Rvar output.}

\item{include_allele_freq}{If TRUE, include allele frequency in .Rvar output.}

\item{show_progress_bar}{If TRUE, print progress bar.}

\item{verbose}{If TRUE, print detail messages.}

\item{logfile}{The log filename. If NULL, print log info on screen.}
}
\value{
a data frame of region_metatable, with region definitions and
filenames of LD matrices and variant information.
}
\description{
Converts PLINK genotype data to LD matrices and SNP info files,
saves LD matrices as .RDS files and SNP info as .Rvar files
}
