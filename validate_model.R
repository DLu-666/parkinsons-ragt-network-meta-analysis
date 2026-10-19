suppressPackageStartupMessages(library(rjags))
suppressPackageStartupMessages(library(metafor))
source('nma_functions.R')
W<-'.'
a<-read.csv(file.path(W,'primary_endpoints.csv'),fileEncoding='UTF-8-BOM')
groups<-split(a,interaction(a$study,a$domain,drop=TRUE))
delta<-numeric();blocks3<-list()
for(z in groups) {
  zz<-smd_study(z)
  if(nrow(z)==2) {
    ref<-escalc(measure='SMD',m1i=z$mean[2]*z$direction[2],m2i=z$mean[1]*z$direction[1],
       sd1i=z$sd[2],sd2i=z$sd[1],n1i=z$n[2],n2i=z$n[1],vtype='UB')
    delta<-c(delta,abs(zz$y-ref$yi),abs(zz$V-ref$vi))
  } else {
    # Rebase from arm1 to arm2; new ordering2,1,3. Linear map must preserve mean and covariance.
    B<-rbind(c(-1,0),c(-1,1));zz2<-smd_study(z,ref=2)
    delta<-c(delta,abs(B%*%zz$y-zz2$y),abs(B%*%zz$V%*%t(B)-zz2$V))
    blocks3[[length(blocks3)+1]]<-zz
  }
  stopifnot(min(eigen(zz$V,symmetric=TRUE)$values)>0)
}
stopifnot(max(delta)<1e-10)
cat('MAX EFFECT/COVARIANCE ERROR:',max(delta),'\n')
# Analytical GLS and MCMC fixed-effect estimates must agree within Monte Carlo error.
zz<-a[a$domain=='6MWT',];bl<-lapply(split(zz,zz$study),smd_study)
spec<-make_model(bl,fixed=TRUE);p<-length(spec$nodes)-1
precision<-diag(.01,p);rhs<-rep(0,p)
for(z in bl){
  X<-matrix(0,length(z$y),p)
  for(j in seq_along(z$y)){
    t1<-match(z$treat[j+1],spec$nodes)-1;t0<-match(z$treat[1],spec$nodes)-1
    if(t1>0)X[j,t1]<-1;if(t0>0)X[j,t0]<- -1
  }
  inv<-solve(z$V);precision<-precision+t(X)%*%inv%*%X;rhs<-rhs+drop(t(X)%*%inv%*%z$y)
}
truth<-drop(solve(precision,rhs))
f<-fit_nma(bl,fixed=TRUE,n_iter=15000,burn=3000,adapt=1000)
err<-colMeans(f$draws)[-1]-truth
mcse<-sqrt(diag(solve(precision)))/sqrt(f$diagnostics$ESS)
stopifnot(all(abs(err)<pmax(6*mcse,.01)))
write.csv(data.frame(parameter=spec$nodes[-1],GLS=truth,MCMC=colMeans(f$draws)[-1],difference=err,MCSE=mcse),file.path(W,'validation_gls.csv'),row.names=FALSE)
saveRDS(list(max_covariance_error=max(delta),fixed_validation=f$diagnostics),file.path(W,'validation.rds'))
cat('VALIDATION_COMPLETE\n');q(status=0)
