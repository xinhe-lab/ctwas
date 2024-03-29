% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/susie_univariate_regression.R
\name{univariate_regression}
\alias{univariate_regression}
\title{Perform Univariate Linear Regression Separately for Columns of X}
\usage{
univariate_regression(
  X,
  y,
  Z = NULL,
  center = TRUE,
  scale = FALSE,
  return_residuals = FALSE
)
}
\arguments{
\item{X}{n by p matrix of regressors.}

\item{y}{n-vector of response variables.}

\item{Z}{Optional n by k matrix of covariates to be included in all
regresions. If Z is not \code{NULL}, the linear effects of
covariates are removed from y first, and the resulting residuals
are used in place of y.}

\item{center}{If \code{center = TRUE}, center X, y and Z.}

\item{scale}{If \code{scale = TRUE}, scale X, y and Z.}

\item{return_residuals}{Whether or not to output the residuals if Z
is not \code{NULL}.}
}
\value{
A list with two vectors containing the least-squares
  estimates of the coefficients (\code{betahat}) and their standard
  errors (\code{sebetahat}). Optionally, and only when a matrix of
  covariates \code{Z} is provided, a third vector \code{residuals}
  containing the residuals is returned.
}
\description{
This function performs the univariate linear
regression y ~ x separately for each column x of X. Each regression
is implemented using \code{.lm.fit()}. The estimated effect size
and stardard error for each variable are outputted.
}
\examples{
set.seed(1)
n = 1000
p = 1000
beta = rep(0,p)
beta[1:4] = 1
X = matrix(rnorm(n*p),nrow = n,ncol = p)
X = scale(X,center = TRUE,scale = TRUE)
y = drop(X \%*\% beta + rnorm(n))
res = univariate_regression(X,y)
plot(res$betahat/res$sebetahat)

}
