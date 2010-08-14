`DoISVA` <-
function(data.m,pheno.v,cf.m,factor.log,pvthCF=0.01,th=0.05,ncomp=NULL){

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

