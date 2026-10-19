suppressPackageStartupMessages(library(rjags))
source('nma_functions.R')
W <- '.'
a <- read.csv(file.path(W,'primary_endpoints.csv'),fileEncoding='UTF-8-BOM')
connected <- function(blocks,pair) {
  reached <- pair[1]
  repeat {old<-reached;for(b in blocks)if(any(b$treat %in% reached))reached<-union(reached,b$treat)
    if(length(old)==length(reached))break}
  pair[2] %in% reached
}
checks<-list();individual<-list();pred<-list();diagall<-list()
for(domain in sort(unique(a$domain))) {
 z<-a[a$domain==domain,];bl<-lapply(split(z,z$study),smd_study)
 f<-readRDS(file.path(W,'fits',paste0('main__',domain,'.rds')))
 # Posterior predictive intervals for a new study contrast, conditional on consistency.
 set.seed(20261013+match(domain,sort(unique(a$domain))))
 tt<-as.matrix(f$samples)[,'tau']
 for(pair in combn(f$nodes,2,simplify=FALSE)) {
  id<-paste(domain,paste(pair,collapse='_'),sep='__')
  direct<-Filter(function(b)all(pair %in% b$treat),bl)
  indirect<-Filter(function(b)!all(pair %in% b$treat),bl)
  dy<-dv<-numeric(0)
  for(b in direct) {
   cc<-as.numeric(b$treat[-1]==pair[1])-as.numeric(b$treat[-1]==pair[2])
   yy<-sum(cc*b$y);vv<-drop(t(cc)%*%b$V%*%cc)
   individual[[paste(id,b$study)]]<-data.frame(domain=domain,study=b$study,label=b$label,treatment1=pair[1],treatment2=pair[2],g=yy,se=sqrt(vv),lower=yy-1.96*sqrt(vv),upper=yy+1.96*sqrt(vv))
   dy<-c(dy,yy);dv<-c(dv,vv)
  }
  np<-f$draws[,pair[1]]-f$draws[,pair[2]]+rnorm(nrow(f$draws),sd=tt)
  pred[[id]]<-data.frame(domain=domain,treatment1=pair[1],treatment2=pair[2],predictive_lower=quantile(np,.025),predictive_upper=quantile(np,.975),direct_studies=length(direct))
  status<-if(length(direct)==0)'Indirect evidence only' else if(!connected(indirect,pair))'No independent indirect path after removing all trials containing the comparison' else 'Node split fitted'
  row<-data.frame(domain=domain,treatment1=pair[1],treatment2=pair[2],direct_studies=length(direct),status=status,direct_mean=NA,direct_lower=NA,direct_upper=NA,indirect_mean=NA,indirect_lower=NA,indirect_upper=NA,difference=NA,difference_lower=NA,difference_upper=NA,tail_probability=NA)
  if(status=='Node split fitted') {
   sp<-make_model(indirect);sp$data$Nd<-length(dy);sp$data$yd<-dy;sp$data$vd<-dv
   code<-sub('\\}\\s*$','dd ~ dnorm(0,0.01)\nfor(h in 1:Nd) { yd[h] ~ dnorm(dd,1/(vd[h]+pow(tau,2))) }\n}',sp$code)
   ini<-lapply(1:4,function(i)list(.RNG.name='base::Wichmann-Hill',.RNG.seed=20261013+i*100L))
   con<-textConnection(code);jm<-jags.model(con,sp$data,inits=ini,n.chains=4,n.adapt=5000,quiet=TRUE);close(con)
   update(jm,10000,progress.bar='none');samp<-coda.samples(jm,c('d','dd','tau'),n.iter=60000,thin=5,progress.bar='none')
   mat<-as.matrix(samp);id1<-match(pair[1],sp$nodes);id2<-match(pair[2],sp$nodes)
   di<-mat[,paste0('d[',id1,']')]-mat[,paste0('d[',id2,']')];dd<-mat[,'dd'];diff<-dd-di
   summ<-function(x)c(mean(x),quantile(x,c(.025,.975)))
   row[1,c('direct_mean','direct_lower','direct_upper')]<-summ(dd)
   row[1,c('indirect_mean','indirect_lower','indirect_upper')]<-summ(di)
   row[1,c('difference','difference_lower','difference_upper')]<-summ(diff)
   row$tail_probability<-2*min(mean(diff>0),mean(diff<0))
   pars<-setdiff(colnames(mat),'d[1]');ess<-coda::effectiveSize(samp[,pars,drop=FALSE]);psrf<-coda::gelman.diag(samp[,pars,drop=FALSE],autoburnin=FALSE,multivariate=FALSE)$psrf[,1]
   diagall[[id]]<-data.frame(id=id,parameter=pars,ESS=ess[pars],PSRF=psrf[pars])
   saveRDS(list(samples=samp,model=code,data=sp$data,summary=row),file.path(W,'fits',paste0('nodesplit__',id,'.rds')))
  }
  checks[[id]]<-row
 }
}
write.csv(do.call(rbind,checks),file.path(W,'node_splitting.csv'),row.names=FALSE)
write.csv(do.call(rbind,individual),file.path(W,'individual_contrasts.csv'),row.names=FALSE)
write.csv(do.call(rbind,pred),file.path(W,'predictive_intervals.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagall),file.path(W,'nodesplit_diagnostics.csv'),row.names=FALSE)
cat('NETWORK_CHECKS_COMPLETE\n')
