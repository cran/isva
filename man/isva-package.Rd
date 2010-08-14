\name{isva}
\alias{isva}
\alias{isva}
\docType{package}
\title{
Independent Surrogate Variable Analysis
}
\description{
Independent Surrogate Variable Analysis is a general algorithm for feature selection in the presence of potentially confounding factors, designed for the analysis of large-scale quantitative genomic data (e.g microarrays). It uses Independent Component Analysis to model the confounding factors and allows for heterogeneity within the phenotype of interest.
}
\details{
\tabular{ll}{
Package: \tab isva\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-08-14\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
Two internal functions perform the dimensionality estimation using approximate Random Matrix Theory (EstDimRMT) and modelling of confounding factors using Independent Component Analysis on the residual variation orthogonal to that of the phenotype of interest (isvaFn). DoISVA allows is the main user function, performing feature selection using independent surrogate variables.
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
