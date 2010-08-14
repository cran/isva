\name{isvaFn}
\alias{isvaFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main engine function for inference of independent surrogate variables (ISVs)}
\description{
This function infers statistically independent surrogate variables by performing Independent Component Analysis on the residual variation orthogonal to that of a phenotype of interest.
}
\usage{
isvaFn(data.m, pheno.v, ncomp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.m}{Data matrix. Rows label features. Columns label samples.}
  \item{pheno.v}{Numeric vector encoding phenotype of interest.}
  \item{ncomp}{Optionally specify number of ISVs to look for. By default will use Approximate Random Matrix Theory to infer this number.}
}
\value{
  A list with following entries:
  \item{n.isv}{Number of inferred ISVs.}
  \item{isv}{Matrix of inferred ISVs.}
}

\references{Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Submitted.}
\author{Andrew E Teschendorff}
\examples{
## see example for DoISVA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

