source('nma_functions.R')
W <- '.'
a <- read.csv(file.path(W,'primary_endpoints.csv'),fileEncoding='UTF-8-BOM')
out <- list()
for(domain in unique(a$domain)) {
 f <- readRDS(file.path('fits',paste0('main__',domain,'.rds')))
 for(z in split(a[a$domain==domain,],a$study[a$domain==domain])) {
  b <- smd_study(z);q <- length(b$y)
  C <- b$V+f$tau[1]^2*(matrix(.5,q,q)+diag(.5,q))
  B <- cbind(-rep(1,q),diag(q));L <- t(B)%*%solve(C,B)
  for(pair in combn(seq_along(b$treat),2,simplify=FALSE)) {
   v <- -L[pair[1],pair[2]];stopifnot(v>0)
   out[[length(out)+1]] <- data.frame(domain=domain,study=b$study,treatment1=b$treat[pair[1]],treatment2=b$treat[pair[2]],conductance=v)
  }
 }
}
write.csv(do.call(rbind,out),file.path(W,'study_conductance.csv'),row.names=FALSE)
cat('Conductance recomputed from saved primary fits; no model refit.\n')
