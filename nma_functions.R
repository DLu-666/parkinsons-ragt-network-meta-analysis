# Local project implementation of the normal-identity contrast model.
# Single pooled within-study SD, exact Hedges correction and unbiased variance.
# Two-arm variances are checked against metafor::escalc(vtype='UB').
# Multi-arm covariance is retained explicitly, not split or approximated by min(SE).
smd_study <- function(a, ref=1L) {
  a <- a[c(ref,setdiff(seq_len(nrow(a)),ref)),]
  n <- a$n; m <- a$mean * a$direction
  nu <- sum(n-1); sp <- sqrt(sum((n-1)*a$sd^2)/nu)
  J <- exp(lgamma(nu/2)-.5*log(nu/2)-lgamma((nu-1)/2))
  g <- J*(m[-1]-m[1])/sp
  A <- matrix(1/n[1],length(g),length(g))+diag(1/n[-1],nrow=length(g))
  kappa <- 1-(nu-2)/(nu*J^2)
  V <- A+kappa*tcrossprod(g)
  list(y=g,V=V,treat=a$treatment,study=a$study[1],label=a$label[1],n=sum(n),df=nu,J=J)
}
make_model <- function(blocks, tau_sd=1, fixed=FALSE, independent=FALSE) {
  K <- length(unique(unlist(lapply(blocks,`[[`,'treat'))))
  nodes <- c('RAGT',sort(setdiff(unique(unlist(lapply(blocks,`[[`,'treat'))),'RAGT')))
  N <- length(blocks); M <- max(lengths(lapply(blocks,`[[`,'y')))
  data <- list(N=N,K=K,m=vapply(blocks,function(z)length(z$y),integer(1)),
    y=matrix(NA,N,M),V=array(NA,c(N,M,M)),R=array(NA,c(N,M,M)),
    tr=matrix(NA_integer_,N,M+1),tau_prec=1/tau_sd^2)
  for(i in seq_along(blocks)) {
    q <- data$m[i];data$y[i,1:q]<-blocks[[i]]$y
    data$V[i,1:q,1:q]<-blocks[[i]]$V
    data$R[i,1:q,1:q]<-matrix(.5,q,q)+diag(.5,q)
    data$tr[i,1:(q+1)]<-match(blocks[[i]]$treat,nodes)
  }
  code <- c('model {',if(fixed)'tau <- 0' else 'tau ~ dnorm(0,tau_prec) T(0,)',
    'd[1] <- 0','for (k in 2:K) { d[k] ~ dnorm(0,0.01) }',
    'for (i in 1:N) {',
    'for (j in 1:m[i]) { mu[i,j] <- d[tr[i,j+1]] - d[tr[i,1]]',
    'for (k in 1:m[i]) { C[i,j,k] <- V[i,j,k] + pow(tau,2) * R[i,j,k] } }',
    'P[i,1:m[i],1:m[i]] <- inverse(C[i,1:m[i],1:m[i]])',
    'y[i,1:m[i]] ~ dmnorm(mu[i,1:m[i]],P[i,1:m[i],1:m[i]])','}','}')
  if(fixed)data$tau_prec<-NULL
  list(code=paste(code,collapse='\n'),data=data,nodes=nodes)
}
fit_nma <- function(blocks,tau_sd=1,fixed=FALSE,seed=20261012,n_iter=30000,burn=10000,adapt=5000) {
  spec<-make_model(blocks,tau_sd,fixed)
  ini<-lapply(1:4,function(i)list(.RNG.name='base::Wichmann-Hill',.RNG.seed=seed+i*100L))
  con<-textConnection(spec$code);on.exit(close(con))
  model<-rjags::jags.model(con,spec$data,inits=ini,n.chains=4,n.adapt=adapt,quiet=TRUE)
  update(model,burn,progress.bar='none')
  samples<-rjags::coda.samples(model,c('d',if(!fixed)'tau'),n.iter=n_iter,thin=5,progress.bar='none')
  draws<-as.matrix(samples);d<-draws[,paste0('d[',seq_along(spec$nodes),']'),drop=FALSE]
  colnames(d)<-spec$nodes
  estimable<-setdiff(colnames(draws),'d[1]')
  ess<-coda::effectiveSize(samples[,estimable,drop=FALSE])
  psrf<-coda::gelman.diag(samples[,estimable,drop=FALSE],autoburnin=FALSE,multivariate=FALSE)$psrf[,1]
  rr<-lapply(combn(spec$nodes,2,simplify=FALSE),function(pair){
    z<-d[,pair[1]]-d[,pair[2]]
    data.frame(treatment1=pair[1],treatment2=pair[2],mean=mean(z),median=median(z),lower=quantile(z,.025),upper=quantile(z,.975),prob_benefit=mean(z>0))
  })
  ranks<-t(apply(d,1,function(x)rank(-x,ties.method='average')))
  ranking<-data.frame(treatment=spec$nodes,SUCRA=(ncol(d)-colMeans(ranks))/(ncol(d)-1),prob_best=colMeans(ranks==1))
  # Marginal normal likelihood deviance, all constants included; DIC at posterior mean.
  tt<-if(fixed)rep(0,nrow(draws)) else draws[,'tau']
  devfun<-function(ds,tau){
    sum(vapply(blocks,function(z){
      q<-length(z$y);C<-z$V+tau^2*(matrix(.5,q,q)+diag(.5,q))
      e<-z$y-(ds[z$treat[-1]]-ds[z$treat[1]])
      q*log(2*pi)+as.numeric(determinant(C,logarithm=TRUE)$modulus)+drop(t(e)%*%solve(C,e))
    },numeric(1)))
  }
  D<-vapply(seq_len(nrow(d)),function(i)devfun(d[i,],tt[i]),numeric(1))
  Dhat<-devfun(colMeans(d),mean(tt));pD<-mean(D)-Dhat
  list(samples=samples,draws=d,nodes=spec$nodes,results=do.call(rbind,rr),ranking=ranking,
       diagnostics=data.frame(parameter=names(ess),ESS=as.numeric(ess),PSRF=as.numeric(psrf)),
       DIC=c(Dbar=mean(D),Dhat=Dhat,pD=pD,DIC=mean(D)+pD),
       tau=if(fixed)c(mean=0,lower=0,upper=0) else c(mean=mean(tt),lower=quantile(tt,.025),upper=quantile(tt,.975)),
       settings=list(seed=seed,chains=4,adapt=adapt,burn=burn,n_iter=n_iter,thin=5,tau_sd=tau_sd,fixed=fixed),model=spec$code)
}
