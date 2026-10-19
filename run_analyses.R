suppressPackageStartupMessages(library(rjags))
source('nma_functions.R')
W<-'.';dir.create(file.path(W,'fits'),showWarnings=FALSE)
a<-read.csv(file.path(W,'all_timepoints.csv'),fileEncoding='UTF-8-BOM')
rob<-read.csv(file.path(W,'ROBUST_RCT_by_outcome.csv'),fileEncoding='UTF-8-BOM')
a$high_risk<-rob$Overall[match(paste(a$study,a$domain),paste(rob$study,rob$domain))]=='High Risk'
stopifnot(!anyNA(a$high_risk))
manifest<-list();allres<-list();alldiag<-list();allrank<-list();allsummary<-list()
jobs<-list(main=list(time='post',filter=function(z)z,tau=1),
  exclude_high_risk=list(time='post',filter=function(z)z[!z$high_risk,],tau=1),
  exclude_median_conversion=list(time='post',filter=function(z)z[z$approximation=='False',],tau=1),
  exclude_combined=list(time='post',filter=function(z)z[z$combination=='False',],tau=1),
  conservative_picelli=list(time='post',filter=function(z)z[!(z$study %in% c(1,2,3)),],tau=1),
  tau_halfnormal05=list(time='post',filter=function(z)z,tau=.5),
  tau_halfnormal15=list(time='post',filter=function(z)z,tau=1.5),
  followup_1month=list(time='1 month',filter=function(z)z,tau=1),
  followup_3months=list(time='3 months',filter=function(z)z,tau=1),
  followup_6months=list(time='6 months',filter=function(z)z,tau=1))
for (jobname in names(jobs)) {
  job<-jobs[[jobname]]
  for(domain in sort(unique(a$domain))) {
    z<-job$filter(a[a$time==job$time & a$domain==domain,])
    id<-paste(jobname,domain,sep='__');ns<-length(unique(z$study));nt<-length(unique(z$treatment))
    status<-if(ns==0)'No eligible data' else if(ns==1)'Single trial only; no pooled NMA' else 'Fitted exploratory NMA'
    manifest[[id]]<-data.frame(analysis=jobname,domain=domain,n_studies=ns,n_nodes=nt,n=sum(z$n),status=status)
    if(ns<2)next
    bl<-lapply(split(z,z$study),smd_study);dest<-file.path(W,'fits',paste0(id,'.rds'))
    if(file.exists(dest)) {f<-readRDS(dest)} else {
      cat('FIT',id,'studies',ns,'\n');flush.console()
      f<-fit_nma(bl,tau_sd=job$tau,seed=20261012+match(domain,sort(unique(a$domain)))*1000)
      saveRDS(f,dest)
    }
    if(max(f$diagnostics$PSRF)>1.01 || min(f$diagnostics$ESS)<1000){
       cat('EXTENDING',id,'\n');f<-fit_nma(bl,tau_sd=job$tau,seed=20261012+match(domain,sort(unique(a$domain)))*1000,n_iter=90000,burn=20000,adapt=10000);saveRDS(f,dest)
    }
    allres[[id]]<-cbind(analysis=jobname,domain=domain,f$results)
    alldiag[[id]]<-cbind(analysis=jobname,domain=domain,f$diagnostics)
    allrank[[id]]<-cbind(analysis=jobname,domain=domain,f$ranking)
    allsummary[[id]]<-data.frame(analysis=jobname,domain=domain,tau_mean=f$tau[1],tau_lower=f$tau[2],tau_upper=f$tau[3],Dbar=f$DIC['Dbar'],pD=f$DIC['pD'],DIC=f$DIC['DIC'])
  }
}
write.csv(do.call(rbind,manifest),file.path(W,'analysis_manifest.csv'),row.names=FALSE)
write.csv(do.call(rbind,allres),file.path(W,'comparative_effects.csv'),row.names=FALSE)
write.csv(do.call(rbind,alldiag),file.path(W,'convergence_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,allrank),file.path(W,'exploratory_rankings.csv'),row.names=FALSE)
write.csv(do.call(rbind,allsummary),file.path(W,'model_fit_summary.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(W,'sessionInfo.txt'))
cat('ALL_ANALYSES_COMPLETE\n');q(status=0)
