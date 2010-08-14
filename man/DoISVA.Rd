\name{DoISVA}
\alias{DoISVA}
\title{Feature selection using independent surrogate variables}
\description{
This function performs feature selection for features associated with a phenotype of interest in the presence of known potential confounding factors subject to uncertainty or measurement error.
}
\usage{
DoISVA(data.m, pheno.v, cf.m, factor.log, pvthCF = 0.01, th = 0.05, ncomp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.m}{Data matrix: rows label features, columns label samples. It is assumed that number of features is much larger than number of samples.}
  \item{pheno.v}{Numeric vector of length equal to number of columns of data matrix. At present categorical phenotypes are not supported.}
  \item{cf.m}{Matrix of confounding factors. Rows label samples, Columns label confounding factors, which may be numeric or categorical.}
  \item{factor.log}{A logical vector of same length as columns of \code{cf.m}. FALSE indicates factor is to be treated as a numeric, TRUE as categorical.}
  \item{pvthCF}{P-value threshold to call a significant association between an independent surrogate variable and a confounding factor. By default this is 0.01.}
  \item{th}{False discovery rate threshold for feature selection. By default this is 0.05.}
  \item{ncomp}{Number of independent surrogate variables to look for. By default this is NULL, and estimation is performed using approximate Random Matrix Theory.}
}
\value{
  A list with following entries:
  \item{lm}{Matrix of feature regression statistics and P-values.}
  \item{qv}{Estimated sorted q-values (False Discovery Rate).}
  \item{spv}{Sorted P-values.}
  \item{rk}{Ranked index of features.}
  \item{isv}{Matrix of selected independent surrogate variables (ISVs).}
  \item{nsv}{Number of selected ISVs.}
  \item{ndeg}{Number of differentially altered features.}
  \item{deg}{Indices of differentially altered features.}
  \item{pvCF}{P-value matrix of associations between factors (phenotype of interest plus confounding factors) and inferred ISVs. Note that this may be a larger set than the selected ISVs.}
  \item{selisv}{Column indices of selected ISVs.}
}
\references{Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Submitted.}
\author{Andrew E Teschendorff}

\examples{

EstDimRMT <- function(data.m){
 ### standardise matrix
 M <- data.m;
 for(c in 1:ncol(M)){
  M[,c] <- (data.m[,c]-mean(data.m[,c]))/sqrt(var(data.m[,c]));
 }
 sigma2 <- var(as.vector(M));
 Q <- nrow(data.m)/ncol(data.m);
 lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
 lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
 delta <- lambdaMAX - lambdaMIN;#  print(delta);
 step <- round(delta/ncol(data.m),3);
 lambda.v <- seq(lambdaMIN,lambdaMAX,by=step);
 dens.v <- vector();
 ii <- 1;
 for(i in lambda.v){
   dens.v[ii] <- (Q/(2*pi*sigma2))*sqrt( (lambdaMAX-i)*(i-lambdaMIN) )/i;
   ii <- ii+1;
 }
 thdens.o <- list(min=lambdaMIN,max=lambdaMAX,step=step,lambda=lambda.v,dens=dens.v);
 C <- 1/nrow(M) * t(M) \%*\% M;
 eigen.o <- eigen(C,symmetric=TRUE);
 estdens.o <- density(eigen.o$values,from=min(eigen.o$values),to=max(eigen.o$values),cut=0);
 intdim <- length(which(eigen.o$values > thdens.o$max));
 return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o));
}

isvaFn <- function(data.m,pheno.v,ncomp=NULL){

  lm.o <- lm(t(data.m) ~ pheno.v);
  res.m <- t(lm.o$res);
  model <- model.matrix(~1+pheno.v);
  if(is.null(ncomp)){
    rmt.o <-  EstDimRMT(res.m)
    ncomp <- rmt.o$dim;
    print(ncomp);
  }
  else {
    print("no need to estimate dimensionality");
  }

  ### perform ICA on residual matrix
  fICA.o <- fastICA(res.m,n.comp=ncomp);
  ### now construct ISV
  tmp.m <- t(fICA.o$A);
  isv.m <- tmp.m;
  for(k in 1:ncol(tmp.m)){
   pv.v <- apply(data.m,1,function(x){summary(lm(tmp.m[,k] ~ x))$coeff[2,4]});
   tmp.s <- sort(pv.v,decreasing=FALSE,index.return=TRUE);
   qv.o <- qvalue(pv.v);
   nsig <- length(which(qv.o$qvalues<0.05));
   if( nsig < 500 ){
     nsig <- 500;
   }
   red.m <- data.m[tmp.s$ix[1:nsig],];
   fICA.o <- fastICA(red.m,n.comp=ncomp);
   cor.v <- abs(cor(tmp.m[,k],t(fICA.o$A)));
   kmax <- which.max(cor.v);
   isv.m[,k] <- t(fICA.o$A)[,kmax];
   print(paste("Done for ISV ",k,sep=""));   
  }
  return(list(n.isv=ncomp,isv=isv.m));
}


## The DoISVA function is currently defined as
DoISVA <- function(data.m,pheno.v,cf.m,factor.log,pvthCF=0.01,th=0.05,ncomp=NULL){

 ### Main ISVA function
 isva.o <- isvaFn(data.m,pheno.v,ncomp);
 ### study pattern of correlation of ISVA components to POI and CFs
 tmp.m <- cbind(pheno.v,cf.m);
 treatfactor <- c(FALSE,factor.log);
 pv.m <- matrix(nrow=ncol(isva.o$isv),ncol=1+ncol(cf.m));
 colnames(pv.m) <- c("POI",colnames(cf.m)); ## POI:phenotype of interest
 for(c in 1:ncol(tmp.m)){
  if(treatfactor[c]==FALSE){
   for(sv in 1:ncol(isva.o$isv)){
    lm.o <- lm(isva.o$isv[,sv] ~ as.numeric(tmp.m[,c]));
    pv.m[sv,c] <- summary(lm.o)$coeff[2,4];   
   }
  }
  else {
   for(sv in 1:ncol(isva.o$isv)){
    lm.o <- lm(isva.o$isv[,sv] ~ as.factor(tmp.m[,c]));
    pv.m[sv,c] <- pf(summary(lm.o)$fstat[1],summary(lm.o)$fstat[2],summary(lm.o)$fstat[3],lower.tail=FALSE);   
   }
  }
 }

 ### selection of ISVs
 selisv.idx <- vector();
 for(sv in 1:nrow(pv.m)){

   ncf <- length(which(pv.m[sv,2:ncol(pv.m)]< pvthCF)) ## pvth=0.01
   minpv <- min(pv.m[sv,2:ncol(pv.m)]);
   phpv <- pv.m[sv,1];
   if(ncf > 0){
     if(minpv < phpv){
       selisv.idx <- c(selisv.idx,sv);
     }
   }
 }
 lm.m <- t(apply(data.m,1,function(x){summary(lm(x ~ pheno.v + isva.o$isv[,selisv.idx]))$coeff[2,3:4]}));
 pv.s <- sort(lm.m[,2],decreasing=FALSE,index.return=TRUE);
 qv.v <- qvalue(pv.s$x)$qvalue;
 ntop <- length(which(qv.v < th));

 if(ntop>0){
  pred.idx <- pv.s$ix[1:ntop];
 }
 else {
  pred.idx <- NULL;
 }
 
 return(list(lm=lm.m,qv=qv.v,spv=pv.s$x,rk=pv.s$ix,isv=isva.o$isv[,selisv.idx],nsv=length(selisv.idx),ndeg=ntop,deg=pred.idx,pvCF=pv.m,selisv=selisv.idx));
 
} ### END OF FUNCTION


### Example
### load in simulated data
data(simdataISVA);
## factors matrix (2 CFs)
factors.m <- cbind(simdataISVA$factors[[1]],simdataISVA$factors[[2]]);
colnames(factors.m) <- c("CF1","CF2");
### Do ISVA
isva.o <- DoISVA(simdataISVA$data,simdataISVA$pheno,factors.m,factor.log=rep(FALSE,2),pvthCF=0.01,th=0.05,ncomp=NULL);
### Evaluation
### modeling of CFs
print(cor(isva.o$isv,factors.m));
### sensitivity
print(length(intersect(isva.o$deg,simdataISVA$deg))/length(simdataISVA$deg));
### PPV (1-false discovery rate)
print(length(intersect(isva.o$deg,simdataISVA$deg))/length(isva.o$deg));



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

