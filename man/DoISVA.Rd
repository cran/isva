\name{DoISVA}
\alias{DoISVA}
\title{Feature selection using independent surrogate variables}
\description{
Given a data matrix and a phenotype of interest, this function performs
feature selection to identify features associated with the phenotype of
interest in the presence of potential confounding factors. The algorithm
first finds the variation in the data matrix not associated with the
phenotype of interest (using a linear model), and subsequently performs
Independent Component Analysis (ICA) on this residual variation
matrix. The number of independent components to be inferred can be
prespecified or estimated using Random Matrix Theory (RMT). Independent
Surrogate Variables (ISVs) are constructed from the independent
components and provide estimates of the effect of confounders on the
data. If potential confounders are unknown (default NULL option) there
will be as many ISVs as there are independent components in the residual
variation space. If potential confounders are known (either exactly or
subject to error/uncertainty) the algorithm will select only those
independent components that correlate with the confounders. If potential
confounders are specified it can happen that ISVA will not select any
ISVs because none of the independent components correlates with the
confounders. In this scenario ISVA should be rerun with the default
(NULL) option. The constructed ISVs are finally included as covariates
in a multivariate regression model to identify features that correlate
with the phenotype of interest independently of the potential
confounders. There are two implementations of ICA offered: JADE and
fastICA. We note that the former will result in the same solution
(therefore deterministic), whereas fastICA may result in convergence to
different local minima for different runs. In the latter case, a
consensus solution is advised if results vary between runs.
}
\usage{
DoISVA(data.m, pheno.v, cf.m = NULL, factor.log, pvthCF = 0.01,
 th = 0.05, ncomp = NULL,icamethod=c("JADE","fastICA"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.m}{Data matrix: rows label features, columns label samples. It is assumed that number of features is much larger than number of samples.}
  \item{pheno.v}{Numeric vector of length equal to number of columns of data matrix. At present only numeric (ordinal) phenotypes are supported, so categorical phenotypes are excluded.}
  \item{cf.m}{Matrix of potential confounding factors. Rows label samples, Columns label confounding factors, which may be numeric or categorical. The default option (NULL) is for the case where potential confounding factors are not known or irrelevant.}
  \item{factor.log}{A logical vector of same length as columns of \code{cf.m}. FALSE indicates factor is to be treated as a numeric, TRUE as categorical.}
  \item{pvthCF}{P-value threshold to call a significant association between an independent surrogate variable and a confounding factor. By default this is 0.01.}
  \item{th}{False discovery rate threshold for feature selection. By default this is 0.05.}
  \item{ncomp}{Number of independent surrogate variables to look for. By
    default this is NULL, and estimation is performed using Random
    MatrixTheory.}
  \item{icamethod}{Method implementing ICA to be used. Must be either
    JADE or fastICA.}  
}
\value{
  A list with following entries:
  \item{spv}{Sorted P-values.}
  \item{rk}{Ranked index of features.}
  \item{qv}{Estimated sorted q-values (False Discovery Rate).}
  \item{ndeg}{Number of differentially altered features.}
  \item{deg}{Indices of differentially altered features.}
  \item{lm}{Matrix of significant feature regression statistics and P-values.}
  \item{isv}{Matrix of selected independent surrogate variables (ISVs).}
  \item{nsv}{Number of selected ISVs.}
  \item{pvCF}{P-value matrix of associations between factors (phenotype of interest plus confounding factors) and inferred ISVs. Note that this may be a larger set than the selected ISVs.}
  \item{selisv}{Column indices of selected ISVs.}
}
\references{Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Teschendorff AE, Zhuang JJ, Widschwendter M. Bioinformatics. 2011 Jun 1;27(11):1496-505.}
\author{Andrew E Teschendorff}

\examples{

### Example

### load in simulated data
data(simdataISVA);
data.m <- simdataISVA$data;
pheno.v <- simdataISVA$pheno;

## factors matrix (two potential confounding factors, e.g chip and cohort)
factors.m <- cbind(simdataISVA$factors[[1]],simdataISVA$factors[[2]]);
colnames(factors.m) <- c("CF1","CF2");

### Estimate number of significant components of variation
rmt.o <- EstDimRMT(data.m);
print(paste("Number of significant components=",rmt.o$dim,sep=""));
### this makes sense since 1 component is associated with the
### the phenotype of interest, while the other two are associated
### with the confounders
ncp <- rmt.o$dim-1 ;

### Do ISVA
### run with the confounders as given
isva.o <- DoISVA(data.m,pheno.v,factors.m,factor.log=rep(FALSE,2),
pvthCF=0.01,th=0.05,ncomp=ncp,icamethod="fastICA");

### Evaluation (ISVs should correlate with confounders)
### modeling of CFs
print(cor(isva.o$isv,factors.m));
### this shows that CFs are reconstructed fairly well

### sensitivity (fraction of detected true positives)
print(length(intersect(isva.o$deg,simdataISVA$deg))/length(simdataISVA$deg));

### PPV (1-false discovery rate)
print(length(intersect(isva.o$deg,simdataISVA$deg))/length(isva.o$deg));

### run not knowing what confounders there are and with ncp=3 say.
isva2.o <- DoISVA(data.m,pheno.v,cf.m=NULL,factor.log=rep(FALSE,2),
pvthCF=0.01,th=0.05,ncomp=3,icamethod="fastICA");

### sensitivity (fraction of detected true positives)
print(length(intersect(isva2.o$deg,simdataISVA$deg))/length(simdataISVA$deg));

### PPV (1-false discovery rate)
print(length(intersect(isva2.o$deg,simdataISVA$deg))/length(isva2.o$deg));


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

