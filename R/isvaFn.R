`isvaFn` <-
function(data.m,pheno.v,ncomp=NULL){

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
  library(qvalue);
  library(fastICA);
  fICA.o <- fastICA(res.m,n.comp=ncomp);

  ### now construct ISV
  tmp.m <- t(fICA.o$A);
  isv.m <- tmp.m;
  for(k in 1:ncol(tmp.m)){
   lm.o <- summary(lm(t(data.m) ~ tmp.m[,k]));
   pv.v <- unlist(lapply(lm.o,function(x){ x$coeff[2,4]}));
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

