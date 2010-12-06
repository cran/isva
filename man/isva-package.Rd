\name{isva}
\alias{isva}
\alias{isva}
\docType{package}
\title{
Independent Surrogate Variable Analysis
}
\description{
Independent Surrogate Variable Analysis is a general algorithm for feature selection in the presence of potential confounding factors, designed for the analysis of large-scale quantitative genomic data (e.g microarrays). It uses Independent Component Analysis (ICA) to model the confounding factors as independent surrogate variables (ISVs). These ISVs are included as covariates in a multivariate regression model to identify features that correlate with a phenotype of interest independently of these confounders. The ICA implementation used is that of the fastICA R-package.
}
\details{
\tabular{ll}{
Package: \tab isva\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2010-12-06\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Two internal functions perform the dimensionality estimation using approximate Random Matrix Theory (EstDimRMT) and modelling of confounding factors using Independent Component Analysis on the residual variation orthogonal to that of the phenotype of interest (isvaFn). DoISVA is the main user function, performing feature selection using independent surrogate variables.
}
\author{
Andrew E Teschendorff
Maintainer:<a.teschendorff@ucl.ac.uk>
}
\references{
Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Submitted.
}
\keyword{multivariate}
\seealso{
}
\examples{
}
