`EstDimRMT` <-
function(data.m){
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
 ## theoretical density
 thdens.o <- list(min=lambdaMIN,max=lambdaMAX,step=step,lambda=lambda.v,dens=dens.v);
 C <- 1/nrow(M) * t(M) %*% M;
 eigen.o <- eigen(C,symmetric=TRUE);
 ## empirical density
 estdens.o <- density(eigen.o$values,from=min(eigen.o$values),to=max(eigen.o$values),cut=0);
 intdim <- length(which(eigen.o$values > thdens.o$max));
 return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o));
}

