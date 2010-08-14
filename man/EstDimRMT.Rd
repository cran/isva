\name{EstDimRMT}
\alias{EstDimRMT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimates dimensionality of data set using Approximate Random Matrix Theory}
\description{
Given the data matrix, it estimates the approximate number of significant components of variation by comparing the observed distribution of spectral eigenvalues to the theoretical one under a Gaussian Orthogonal Ensemble (GOE).
}
\usage{
EstDimRMT(data.m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.m}{Data matrix. Rows label features, Columns samples.}
}
\value{
  A list with following objects
  \item{cor}{Data covariance matrix.}
  \item{dim}{Estimated intrinsic dimensionality of data.}
  \item{estdens}{Empirical density of eigenvalues.}
  \item{thdens}{Theoretical density of eigenvalues.}
}
\references{Random matrix approach to cross correlations in financial data. Plerou et al. Physical Review E (2002), Vol.65.}
\author{Andrew E Teschendorff}
\examples{
## see example for DoISVA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

