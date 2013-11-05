\name{isva}
\alias{isva}
\docType{package}
\title{
Independent Surrogate Variable Analysis
}
\description{
Independent Surrogate Variable Analysis is an algorithm for feature selection in the presence of potential confounding factors, specially designed for the analysis of large-scale high-dimensional quantitative genomic data (e.g microarrays). It uses Independent Component Analysis (ICA) to model the confounding factors as independent surrogate variables (ISVs). These ISVs are included as covariates in a multivariate regression model to subsequently identify features that correlate with a phenotype of interest independently of these confounders. The ICA implementation used is that of the fastICA R-package.
}
\details{
\tabular{ll}{
Package: \tab isva\cr
Type: \tab Package\cr
Version: \tab 1.7\cr
Date: \tab 2012-09-25\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
There are two internal functions. One function (EstDimRMT) performs the dimensionality estimation using a Random Matrix Theory approximation. The other function (isvaFn) is the main engine function and performs the modelling of confounding factors using Independent Component Analysis (ICA). Briefly, ICA is applied on the residual variation orthogonal to that of the phenotype of interest. DoISVA is the main user function, performing feature selection using the constructed independent surrogate variables as covariates.
}
\author{
Andrew E Teschendorff
Maintainer:<a.teschendorff@ucl.ac.uk>
}
\references{
Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Teschendorff AE, Zhuang JJ, Widschwendter M. Bioinformatics. 2011 Jun 1;27(11):1496-505.
}
\keyword{multivariate}

