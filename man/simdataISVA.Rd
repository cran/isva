\name{simdataISVA}
\alias{simdataISVA}
\docType{data}
\title{Simulated data for ISVA}
\description{
  A synthetic data set of 2000 features and 50 samples with a binary phenotype
  and two confounding factors. Relative effect size of confounding
  factors (CFs) to that of phenotype of interest is 4. For further details please
 see reference.
}
\usage{simdataISVA}
\format{
  This synthetic data set is a list object containing the following
  elements: (i) data is the data matrix (2000 features, 50 samples),
  (ii) pheno is a binary phenotype vector, (iii) factors is a list of
  length two containing the two binary confounding factors, (iv) deg is
  the index vector of those truly differentially "expressed" features, (v)
  degL is a list of index vectors for features truly differentially
  altered (first element,degL[[1]]=deg) and those features affected by CFs (2nd and
  3rd elements).
}
\references{Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Teschendorff AE, Zhuang JJ, Widschwendter M. Bioinformatics. 2011 Jun 1;27(11):1496-505.}
\keyword{datasets}