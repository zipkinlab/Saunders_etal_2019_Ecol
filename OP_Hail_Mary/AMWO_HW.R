# Run JAGS
setwd("/mnt/scratch/rossmans/AMWO_HM")
laod(AMWO_HM.RData)
AMWO.combo.jags<-NA
class(AMWO.combo.jags)<-"try-error"
breaks<-0
while(class(run.it)=="try-error"){
  #if (breaks>20){break}
  inits <- function(){list(H=H.inits,pi.sex=pi.sex.inits, sa.x=sa.x.inits, ss.x=ss.x.inits, n=n.inits, f.sd=f.sd.inits,f.x=f.x.inits, F.x=F.x.inits)} #,report=report.inits)}                 
  breaks<-breaks+1
  AMWO.combo.jags <- try(jagsUI(bugs.data, inits=NULL, parameters, "Lincoln.Brownie.30March.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE)) #changed from inits=NULL
}
save.image(paste("SUCCESS",floor(runif(1,0,10000)),".RData",sep=""))
