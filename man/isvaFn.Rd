\name{isvaFn}
\alias{isvaFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main engine function for inference of independent surrogate variables (ISVs)}
\description{
This is the main engine function which infers the statistically
independent surrogate variables (ISVs) by performing Independent
Component Analysis (ICA) on the residual variation matrix. It uses
either the ICA implementation of JADE or the one from the fastICA R-package. The residual variation matrix reflects the variation orthogonal to that of a phenotype of interest and is inferred using a linear model.
}
\usage{
isvaFn(data.m, pheno.v, ncomp = NULL,icamethod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.m}{Data matrix. Rows label features. Columns label samples.}
  \item{pheno.v}{Numeric vector encoding phenotype of interest.}
  \item{ncomp}{Optionally specify number of ISVs to look for. By default
    will use Approximate Random Matrix Theory to infer this number.}
  \item{icamethod}{The ICA method to be used. Input value is taken from DoISVA.}  
}
\value{
  A list with following entries:
  \item{n.isv}{Number of inferred ISVs.}
  \item{isv}{Matrix of inferred ISVs.}
}

\references{Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Teschendorff AE, Zhuang JJ, Widschwendter M. Bioinformatics. 2011 Jun 1;27(11):1496-505.}
\author{Andrew E Teschendorff}
\examples{
## see example for DoISVA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

