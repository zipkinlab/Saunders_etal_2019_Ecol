#--------------------------------------------------------------------
#American woodcock integrated population model
#Code by Saunders et al.
#January 2017 - December 2017
#Submitted & in review at Ecology
#--------------------------------------------------------------------

#load jagsUI
library(jagsUI)

# Load data objects for HPCC run
load(file="AMWO_harvest_pi_Oct2017.Rda")
#load m-array for IPM and release array
load(file="marrayAMWO.Rda")
load(file="relAMWO.Rda")

n.yrs = 50

sgs <- read.csv("SGS-indices.csv", header=TRUE)
region <- sgs[,1]
region <- factor(region, levels=rev(levels(region)))
region <- as.numeric(region)

sgs <- sgs[,-c(1:2)]

#Add 5 years of NA's to match the IPM output years
sgs.matrix <- matrix(NA,nrow=25,ncol=54)
sgs.matrix[,6:54] <- as.matrix(sgs)
sgs <- t(sgs.matrix)

#------------#
#-BUGS Model-#
#------------#

sink("AMWO.IPM.Dec2017.jags")
cat("
    model {
    
    # Priors and constraints for population means and variances
    #first dimension is t, second dimension is season, third dimension is age-sex class [1=juv,2=male,3=female], last is population [1=eastern,2=central]
    #NOTE: Values below are SCALED--divided by 10000 for all    
    
    # prior for initial pop sizes in spring               # informed prior based on Lincoln estimates
    n[1,1,3,1] ~ dunif(100, 800)                          # initial size for female eastern
    n[1,1,3,2] ~ dunif(100, 800)                          # initial size for female central
    n[1,1,2,1] ~ dunif(100, 800)                          # initial size for male eastern
    n[1,1,2,2] ~ dunif(100, 800)                          # initial size for male central
    n[1,1,1,1] <- n[1,1,3,1] * F[1,1]                     # initial size for juv eastern is function of female initial size
    n[1,1,1,2] <- n[1,1,3,2] * F[1,2]                     # initial size for juv central is function of female initial size
    
    # round initial sizes
    N[1,1,3,1] <- round(n[1,1,3,1])
    N[1,1,3,2] <- round(n[1,1,3,2])
    N[1,1,2,1] <- round(n[1,1,2,1])
    N[1,1,2,2] <- round(n[1,1,2,2])
    N[1,1,1,1] <- round(n[1,1,1,1])
    N[1,1,1,2] <- round(n[1,1,1,2])
    
    for (p in 1:2){                                             # 1 Eastern, 2 Central
    # priors for fecundities
    F.x[p] ~ dunif(0, 4)                                        # logical bounds on fecundity (F) with 4 egg clutch
    F.sd[p] ~ dunif(0.01, 1)                                    # Fecundity cannot be <0 or >4
    F.tau[p] <- pow(F.sd[p], -2)
    
    for (t in 1:yrs){                                           #yrs = 53
    F[t,p] ~ dnorm(F.x[p], F.tau[p])T(0,)                       
    
    # deriving pi.sex for use in balance equations
    pi.sex[t,1,p] <- pi[t,2,p]/(pi[t,2,p]+pi[t,3,p])            # male to female proportion (pulling pi's from harvest model run)
    pi.sex[t,2,p] <- 1-pi.sex[t,1,p]                            # female to male proportion
    
    N.tot.sum[t,p] <- N[t,2,1,p] + N[t,2,2,p] + N[t,2,3,p]          #deriving total population size for each pop post summer
    N.tot.sp[t,p] <- N[t,1,1,p] + N[t,1,2,p] + N[t,1,3,p]           #deriving total population size for each pop spring
    
    } #t
    
    for (c in 1:3){                                             #Priors for survival and recovery rates for each age-sex class and population
    sa.x[c,p] ~ dunif(0,1)                                      #annual survival
    sa.mu[c,p] <- logit(sa.x[c,p])        
    ss.x[c,p] ~ dunif(0,1)                                      #summer survival
    ss.mu[c,p] <- logit(ss.x[c,p])
    f.x[c,p] ~ dunif(0,0.1)                                     #recovery rate
    f.mu[c,p] <- logit(f.x[c,p])               
    
    sa.sd[c,p] ~ dunif(0.05,2)                                  # Priors for SDs of survival and recovery rates
    ss.sd[c,p] ~ dunif(0.05,2) 
    f.sd[c,p] ~ dunif(0.05,2) 
    sa.tau[c,p] <- pow(sa.sd[c,p],-2)                           # Express as precision
    ss.tau[c,p] <- pow(ss.sd[c,p],-2) 
    f.tau[c,p] <- pow(f.sd[c,p],-2) 
    
    ## Generate annual parameter estimates
    for (t in 1:yrs){
    eps.sa[t,c,p] ~ dnorm(0,sa.tau[c,p])
    eps.ss[t,c,p] ~ dnorm(0,ss.tau[c,p])
    eps.f[t,c,p] ~ dnorm(0,f.tau[c,p])
    
    logit(sa[t,c,p]) <- sa.mu[c,p] + eps.sa[t,c,p]             # annual survival (by year, age-sex class, pop)
    logit(ss[t,c,p]) <- ss.mu[c,p] + eps.ss[t,c,p]             # seasonal survival (summer)
    logit(f[t,c,p]) <- f.mu[c,p] + eps.f[t,c,p]                # Brownie recovery rate
    
    } # close t
    
    ## Balance equations, assuming estimating N in both seasons
    # Initial pop sizes in summer [t, s, c, p]
    N[1,2,c,p] <- trunc(N[1,1,c,p]*ss[1,c,p])                  # pop size in late summer is function of summer survival 
    } #c
    
    for (t in 2:yrs){
    N[t,1,1,p] <- trunc(N[t,1,3,p] * F[t,p])                   # pop size of juv in spring = adF * mean fecundity
    N[t,2,1,p] <- trunc(N[t,1,1,p] * ss[t,1,p])                # pop size of juv in late summer = juv * summer survival
    sw[t,1,p] <- sa[t-1,1,p]/ss[t,1,p]                         # derived winter survival for juveniles (annual/summer)
    
    for (c in 2:3){
    sw[t,c,p] <- sa[t-1,c,p]/ss[t,c,p]                                                        # derived winter survival for adults
    N[t,1,c,p] <- trunc(N[t-1,2,c,p]*sw[t,c,p] + pi.sex[t,c-1,p]*N[t-1,2,1,p]*sw[t,1,p])      # pop size of adults in spring
    N[t,2,c,p] <- trunc(N[t,1,c,p]*ss[t,c,p])                                                # pop size of adults in late summer 
    
    } #c
    } #t 
    
    # Observation process, recoveries and harvest data
    # Recoveries
    # Note: releases MUST be provided as data; saved as relAMWO
    for (c in 1:3){
    for (t in 1:yrs){
    marrayAMWO[t,1:(yrs+1),1,c,p] ~ dmulti(pr[t,,1,c,p], rel[t,1,c,p])        # Apr-Jun releases     
    marrayAMWO[t,1:(yrs+1),2,c,p] ~ dmulti(pr[t,,2,c,p], rel[t,2,c,p])        # Jul-Sep releases
    
    ## Define cell probabilities of m-arrays
    # Main diagonal--direct recovery in season [t]
    pr[t,t,1,c,p] <- ss[t,c,p] * f[t,c,p]                  # spring banded birds must survive to start of hunting season
    pr[t,t,2,c,p] <- f[t,c,p]                              # survival assumed 1 if banded Jul-Sep
    
    # monitor cumulative survival to start of next hunting season (used in subsequent diagonals)
    cumS[t,t,1,c,p] <- ss[t,c,p] * sa[t,c,p]
    cumS[t,t,2,c,p] <- sa[t,c,p]
    } # t
    } # c
    
    for (t in 1:yrs){
    # Above main diagonal--indirect recovery in season [k > t]
    # All birds are adults now, no differences between Apr-Jun and Jul-Sep either
    for (k in (t+1):yrs){                                                                           # k loop represents next year (above main diagonal)
    # recoveries
    pr[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi.sex[t,1,p]*f[k,2,p] + pi.sex[t,2,p]*f[k,3,p])          # juvs as mixture of AdM & AdF
    pr[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi.sex[t,1,p]*f[k,2,p] + pi.sex[t,2,p]*f[k,3,p])          # pi.sex[1] is proportion male; pi.sex[2] is proportion female as esimated from harvest model run      
    pr[t,k,1,2,p] <- cumS[t,k-1,1,2,p] * f[k,2,p] 
    pr[t,k,2,2,p] <- cumS[t,k-1,2,2,p] * f[k,2,p] 
    pr[t,k,1,3,p] <- cumS[t,k-1,1,3,p] * f[k,3,p] 
    pr[t,k,2,3,p] <- cumS[t,k-1,2,3,p] * f[k,3,p] 
    # monitor cumulative survival to start of next hunting period
    cumS[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi.sex[t,1,p]*sa[k,2,p] + pi.sex[t,2,p]*sa[k,3,p])      # juvs as mixture AdM & AdF
    cumS[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi.sex[t,1,p]*sa[k,2,p] + pi.sex[t,2,p]*sa[k,3,p]) 
    cumS[t,k,1,2,p] <- cumS[t,k-1,1,2,p] * sa[k,2,p]
    cumS[t,k,2,2,p] <- cumS[t,k-1,2,2,p] * sa[k,2,p]
    cumS[t,k,1,3,p] <- cumS[t,k-1,1,3,p] * sa[k,3,p]
    cumS[t,k,2,3,p] <- cumS[t,k-1,2,3,p] * sa[k,3,p]
    } #k
    
    # Left of main diag
    for (l in 1:(t-1)){                               # l loop for previous year (left of main diagonal)                    
    pr[t,l,1,1,p] <- 0
    pr[t,l,1,2,p] <- 0
    pr[t,l,1,3,p] <- 0
    pr[t,l,2,1,p] <- 0
    pr[t,l,2,2,p] <- 0
    pr[t,l,2,3,p] <- 0
    } #l
    
    # Last column: probability of non-recovery
    pr[t,(yrs+1),1,1,p] <- 1-sum(pr[t,1:yrs,1,1,p])
    pr[t,(yrs+1),1,2,p] <- 1-sum(pr[t,1:yrs,1,2,p])
    pr[t,(yrs+1),1,3,p] <- 1-sum(pr[t,1:yrs,1,3,p])
    pr[t,(yrs+1),2,1,p] <- 1-sum(pr[t,1:yrs,2,1,p])
    pr[t,(yrs+1),2,2,p] <- 1-sum(pr[t,1:yrs,2,2,p])
    pr[t,(yrs+1),2,3,p] <- 1-sum(pr[t,1:yrs,2,3,p])
    } #t
    
    # Bringing in harvest proportion data (H) estimated from harvest model separately
    for (t in 1:yrs){ 
    for (c in 1:3){
    H[t,c,p] ~ dbin(h[t,c,p], N[t,2,c,p])
    
    h[t,c,p] <- f[t,c,p]/report[t]                    # harvest rate (h) is recovery rate (f) divided by reporting rate (report)                                                      
    
    } #c                                              
    } #t           
    } #p
    
    # Incorporating SGS data
    for (p in 1:2){
    
    tau.sgs[p] <- pow(sd.sgs[p], -2)
    sd.sgs[p] ~ dunif(1, 10)
    
    mu[1,p] ~ dnorm(0,0.37)T(0,)
    r[1,p] ~ dnorm(1, 0.37)T(0,)
    
    } #p
    
    for (t in 1:(length(sgs[,1])-1)){
    for (s in 1:length(sgs[1,])){
    
    sgs[t,s] ~ dnorm(mu[t,region[s]], tau.sgs[region[s]])
    } #s
    } #t
    
    for (t in 2:(length(sgs[,1])-1)){
    for (p in 1:2){
    
    mu[t,p] <- r[t,p] * mu[t-1,p]
    r[t,p] <- N[t,1,2,p]/N[t-1,1,2,p]
    
    } #p
    } #t
    
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data
# Data to be pulled in from separate harvest model
H <- array(NA, dim = c(53,3,2))                   # Harvest proportion data
pi <- array(NA, dim=c(53,3,2))                    # male:female proportion data

for(c in 1:3){
  for(p in 1:2){
    for(t in 2:(n.yrs+1)){
      pi[t,c,p] <- AMWO.state.harvest.pi$mean$pi[t-1,c,p]            #pull in pi estimates
      H[t,c,p] <- round(AMWO.state.harvest.pi$mean$H[t-1,c,p])       #pull in harvest estimates
    }
  }
}

#replace missing H and pi values for last 2 years and first year with means of last 5 years
H[1,1,] <- round(mean(H[2:5,1,], na.rm = TRUE))
H[1,2,] <- round(mean(H[2:5,2,], na.rm = TRUE))
H[1,3,] <- round(mean(H[2:5,3,], na.rm = TRUE))
H[52:53,1,] <- round(mean(H[47:51,1,], na.rm = TRUE))
H[52:53,2,] <- round(mean(H[47:51,2,], na.rm = TRUE))
H[52:53,3,] <- round(mean(H[47:51,3,], na.rm = TRUE))
pi[1,1,] <- mean(pi[2:5,1,], na.rm = TRUE)
pi[1,2,] <- mean(pi[2:5,2,], na.rm = TRUE)
pi[1,3,] <- mean(pi[2:5,3,], na.rm = TRUE)
pi[52:53,1,] <- mean(pi[47:51,1,], na.rm = TRUE)
pi[52:53,2,] <- mean(pi[47:51,2,], na.rm = TRUE)
pi[52:53,3,] <- mean(pi[47:51,3,], na.rm = TRUE)

#scale harvest estimates (divide by 10000) from harvest model output
H <- round(H/10000)

#reporting rate vector. Toll free bands start in 2001 (using reporting rate of 0.8 as in literature, 0.506 for all years prior to that)
report<-c(0.506,rep(0.506,37),0.8,rep(0.8,14))

bugs.data <- list(yrs=dim(marrayAMWO)[1], marrayAMWO=marrayAMWO, rel=relAMWO, H=H, pi=pi, report=report, sgs=sgs, region=region)  

# inits
n.inits <- array(NA, dim=c(1,1,3,2))     #scaled these inits by dividing by 10000
n.inits[1,1,3,1] <- 300
n.inits[1,1,3,2] <- 400
n.inits[1,1,2,1] <- 250
n.inits[1,1,2,2] <- 300

inits <- function(){list(
  n=n.inits
)}               

# Parameters monitored
parameters <- c("sa", "f", "ss", "F","N","sa.x", "ss.x", "f.x", "F.x", "sa.sd", "ss.sd", "f.sd", "F.sd", "N.tot.sp", "N.tot.sum", "sw", "r", "mu", "sd.N", "sd.sgs")  # include h?

# MCMC settings (all params converge with these settings upon quick glance)
ni <- 400000
nt <- 10
nb <- 350000
nc <- 3
na <- 10000

# Run JAGS
AMWO.IPM.jags <- jagsUI(bugs.data, inits=inits, parameters, "AMWO.IPM.Dec2017.jags", n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, store.data=TRUE) #, parallel=TRUE)
