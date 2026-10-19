W<-'.'
domains<-c('motor','balance','mobility','fog','qol','10MWT','6MWT')
for (di in seq_along(domains)) {
 f<-readRDS(file.path(W,'fits',paste0('main__',domains[di],'.rds')))
 mats<-lapply(f$samples,as.matrix);pars<-setdiff(colnames(mats[[1]]),'d[1]')
 pdf(file.path(W,paste0('trace_',domains[di],'.pdf')),width=8,height=10.5,family='Times')
 par(mfrow=c(length(pars),2),mar=c(2.3,3.1,1.6,.5),oma=c(1,0,2,0))
 for(pa in pars) {
  vals<-lapply(mats,function(m)m[,pa]);limits<-range(unlist(vals))
  plot(vals[[1]],type='l',col='#005a80',ylim=limits,xlab='',ylab=pa,main='Retained chain draws',cex.main=.85)
  for(k in 2:4)lines(vals[[k]],col=c('#005a80','#bd6d16','#5d8232','#8d467d')[k])
  den<-lapply(vals,density);plot(den[[1]],col='#005a80',main='Posterior density',xlab='',ylab='',xlim=range(vapply(den,function(d)range(d$x),numeric(2))),ylim=c(0,max(vapply(den,function(d)max(d$y),numeric(1)))),cex.main=.85)
  for(k in 2:4)lines(den[[k]],col=c('#005a80','#bd6d16','#5d8232','#8d467d')[k])
 }
 mtext(paste0('Supplementary Figure 2',letters[di],'. ',domains[di],': main model diagnostics'),outer=TRUE,cex=1.1)
 dev.off()
}
