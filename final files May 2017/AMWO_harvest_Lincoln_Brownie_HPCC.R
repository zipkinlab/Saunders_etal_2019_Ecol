#--------------------------------------------------------------------------------------------------------------
# Estimating annual harvest of AMWO during duck stamp years using annual sales as a covariate
#--------------------------------------------------------------------------------------------------------------

#rm(list=ls(all=TRUE)) # clear memory

# load R2WinBUGS, set working directory for data and output, and identify WinBUGS location for your computer
library(jagsUI)

# LOAD DATA #
#harvest <- read.csv("AMWO harvest 2.csv",header=TRUE)
#dim(harvest)
#harvest[1,]

# read data into 51x2 matrices for jags (col 1 = Eastern, 2 = Central)
# row 1 is 2013, row 51 is 1963 (reverse time)
#HIP <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)
#DSS <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)
#stamps <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)

#for (i in 1:51){
  #HIP[i,2] <- harvest[i,13] 
  #HIP[i,1] <- harvest[(i+51),13] 
  #DSS[i,2] <- harvest[i,12] 
  #DSS[i,1] <- harvest[(i+51),12]}

#head(HIP)

##############################
# scale harvest estimates by dividing by 10000
##############################
#HIP <- HIP/10000
#DSS <- DSS/10000

# use mean and sd from first 10 years HIP (2004-2013) for year 1 priors. NOTE: Scaled now
#mean(HIP[1:10,1])  #eastern
#1/(sd(HIP[1:10,1]))^2 # precision
#mean(HIP[1:10,2])  #central
#1/(sd(HIP[1:10,2]))^2

# use max duck stamps as lower bound on max hunters
#max(stamps[,1],na.rm=TRUE)  #eastern
#max(stamps[,2],na.rm=TRUE)  #central

#save(HIP, file="HIP.Rda")
#save(DSS, file="DSS.Rda")

########################################################
# Load WING and HARVEST data #
########################################################
#harvest <- read.csv("AMWO harvest.csv",header=TRUE)

#bring in year
#clean<-matrix(NA,nrow=102,ncol=1)
#clean<-data.frame(clean)
#clean[,1]<- harvest$Season

#bring in population
#clean[,2]<- harvest$Pop

#bring in age classes of harvest
#clean[,3]<-harvest$AF
#clean[,4]<-harvest$AM
#clean[,5]<-harvest$AU
#clean[,6]<-harvest$IF
#clean[,7]<-harvest$IM
#clean[,8]<-harvest$IU

#this is the sum of all age-sex wing classes. 
#clean[,9] <- rowSums(clean[,c(3:4,6:8)])    # NOTE: Adjust if don't want all unknown groups

#colnames(clean)<-c("Year","Pop","AdF","AdM","AdU","JuvF","JuvM","JuvU","Total")

# Now need to make wings.age and wings.sex datasets
#clean[,10] <- rowSums(clean[,3:5])   # adult wings column; adjust categories included if needed
#clean[,5] <- rowSums(clean[,6:8])   # juvenile wings column; adjust categories included if needed

#colnames(clean)<-c("Year","Pop","AdF","AdM","JuvTot","JuvF","JuvM","JuvU","Total")

# wings<-array(NA, dim=c(51,4,2))                # dims: t, ages (1=juv; 2=adult male; 3=adult female; 4=total wings), p (1=eastern; 2=central)
# wings[,1,1]<-clean$JuvTot[clean$Pop=="E"]      #juvenile eastern
# wings[,2,1]<-clean$AdM[clean$Pop=="E"]         #adult male eastern
# wings[,3,1]<-clean$AdF[clean$Pop=="E"]         #adult female eastern
# wings[,4,1]<-clean$Total[clean$Pop=="E"]       # total wings eastern
# wings[,1,2]<-clean$JuvTot[clean$Pop=="C"]      #juvenile central
# wings[,2,2]<-clean$AdM[clean$Pop=="C"]         # adult male central
# wings[,3,2]<-clean$AdF[clean$Pop=="C"]         # adult female central
# wings[,4,2]<-clean$Total[clean$Pop=="C"]       # total wings central

# replace NAs in 1997 and 1998
# wings[35:36,1,1] <- round(mean(wings[,1,1], na.rm=TRUE))
# wings[35:36,2,1] <- round(mean(wings[,2,1], na.rm=TRUE))
# wings[35:36,3,1] <- round(mean(wings[,3,1], na.rm=TRUE))
# wings[35:36,4,1] <- round(mean(wings[,4,1], na.rm=TRUE))
# wings[35:36,1,2] <- round(mean(wings[,1,2], na.rm=TRUE))
# wings[35:36,2,2] <- round(mean(wings[,2,2], na.rm=TRUE))
# wings[35:36,3,2] <- round(mean(wings[,3,2], na.rm=TRUE))
# wings[35:36,4,2] <- round(mean(wings[,4,2], na.rm=TRUE))

#save(wings, file="wings.Rda")

# Load data objects for HPCC run
#load(file="HIP.Rda")
#load(file="DSS.Rda")
#load(file="wings.Rda")

#############################################################
# specify jags model to predict total annual harvest
#############################################################
#sink("AMWO.harvest.wings.jags")
#cat("
    #model {
    
    # [,1] is eastern pop, [,2] is central population
    # prior for year 1 Harvest (reverse time, [1,]=2013)
    #harvest[1] ~ dnorm(16,9.551e-2) # mean, 1/sd2 from 2004-2013 data, eastern SCALED
    #harvest[2] ~ dnorm(44,1.822e-2) # central SCALED
    
    # round harvest estimates for first years
    #Harvest[1,1] <- round(harvest[1])
    #Harvest[1,2] <- round(harvest[2])
    
    # priors for annual variation in true harvest (sigma.H) 
    # and estimated harvest (sd.HIP, sd.DSS)
    #for (p in 1:2){
    #sigma.H[p] ~ dunif(1,20)  # changed from 200000 for scaling
    #tau.H[p] <- pow(sigma.H[p],-2)
    #sd.HIP[p] ~ dunif(1,20)  # changed from 200000 for scaling
    #prec.HIP[p] <- pow(sd.HIP[p],-2)
    #sd.DSS[p] ~ dunif(1,20)  # changed from 200000 for scaling
    #prec.DSS[p] <- pow(sd.DSS[p],-2)
    #} #p
    
    # State process (estimate true Harvest in years 2-51)
    #for (p in 1:2){
    #for (i in 2:yrs){
    #eps[i,p] ~ dnorm(0,tau.H[p])
    #Harvest[i,p] <- trunc(Harvest[i-1,p] + eps[i,p]) # autoregressive model
    #} #i
    #} #p
    
    # observation process, recoveries and harvest data
    #for (p in 1:2){
    #for (i in 1:15){
    #HIP[i,p] ~ dnorm(Harvest[i,p],prec.HIP[p])} #i  end HIP years
    #for (j in 13:50){
    #DSS[j,p] ~ dnorm(Harvest[j,p]/cor[p],prec.DSS[p])} #j end DSS yrs
    #} # end p loop
    
    # Determining age-sex class proportions (pi's)
    #for (p in 1:2){
    #for (t in 1:yrs){
    #pi[t,1,p] ~ dunif(0.1, 0.8)
    #pi[t,2,p] ~ dunif(0.1, 0.8)
    #pi[t,3,p] <- 1-(pi[t,1,p]+pi[t,2,p])
    
    # Summarize known wing samples by age and sex
    # populations (1 eastern, 2 central)
    #wings.class[t,,p] ~ dmulti(pi[t,,p], wings.total[t,p])          # pi is 3 proportions: 1=juvs, 2=males, 3=females    
    #H[t,1,p] <- pi[t,1,p]*Harvest[rev[t],p]
    #H[t,2,p] <- pi[t,2,p]*Harvest[rev[t],p]
    #H[t,3,p] <- pi[t,3,p]*Harvest[rev[t],p]
    #} #t
    #} #p
    #} # end jags model
    #",fill = TRUE)
#sink()

# Bundle data
#wings.total <- wings[,4,]
#wings.class <- wings[,1:3,]

#Reverse time vector
#rev <- seq(51,1,-1)
#cor <- c(3.1, 4.1)
#bugs.data <- list(yrs=51, HIP=HIP, DSS=DSS, wings.total=wings.total, wings.class=wings.class,
                  #rev=rev, cor=cor)

# Initial values (not all priors listed)

#pi.inits <- array(dim=c(51,3,2))
#for (c in 1:3){
  #for (p in 1:2){
    #for (t in 1:51){
      #pi.inits[t,1,p] <- 0.46
      #pi.inits[t,2,p] <- 0.21
      #pi.inits[t,3,p] <- NA
    #}}}

#inits <- function(){list(pi=pi.inits, harvest=c(13,37))}   #harvest inits based on priors (means of 2004-2013), scaled from 130000 and 370000

# Parameters monitored
#parameters <- c("Harvest","sd.DSS","sd.HIP","sigma.H","H","pi")  

# MCMC settings
#ni <- 200000
#nt <- 5          
#nb <- 150000 
#nc <- 3

# call jags
#AMWO.harvest.wings.scaled.HPCC <- jagsUI(bugs.data, inits=inits, parameters, "AMWO.harvest.wings.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
#save(AMWO.harvest.wings.scaled.HPCC, file="AMWO_harvest_wings_scaled_model_HPCC.Rda")
load(file="AMWO_harvest_wings_scaled_model_HPCC.Rda")
#print(AMWO.harvest.wings.scaled.HPCC, digits=3)

load(file="marrayAMWO.Rda")
load(file="relAMWO.Rda")

#------------#
#-BUGS Model-#
#------------#

sink("Lincoln.Brownie.Simple2.6April.jags")
cat("
    model {
    
    # Priors and constraints for population means and variances
    
    # prior for initial pop sizes in spring               # informed prior based on Lincoln estimates
    n[1,1,3,1] ~ dunif(100, 800)                  # initial size for female eastern (dnorm(600000,1E-10)). Scaled--divided by 10000 for all
    n[1,1,3,2] ~ dunif(100, 800)                  # initial size for female central (dnorm(720000,1E-10))
    n[1,1,2,1] ~ dunif(100, 800)                  # initial size for male eastern (dnorm(530000,1E-10))
    n[1,1,2,2] ~ dunif(100, 800)                  # initial size for male central (dnorm(615000,1E-10))
    n[1,1,1,1] <- n[1,1,3,1] * F[1,1]                     # initial size for juv eastern
    n[1,1,1,2] <- n[1,1,3,2] * F[1,2]                     # initial size for juv central
    
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
    F.sd[p] ~ dunif(0.01, 1) 
    F.tau[p] <- pow(F.sd[p], -2)
    
    for (t in 1:yrs){
    F[t,p] ~ dnorm(F.x[p], F.tau[p])T(0,)                       # Fecundity cannot be <0 or >4

    # deriving pi.sex for use in balance equations
    pi.sex[t,1,p] <- pi[t,2,p]/(pi[t,2,p]+pi[t,3,p])            # male to female proportion (pulling pi's from harvest model run above)
    pi.sex[t,2,p] <- 1-pi.sex[t,1,p]                            # female to male proportion
    
    } #t
    
    for (c in 1:3){                       
    sa.x[c,p] ~ dunif(0,1) 
    sa.mu[c,p] <- logit(sa.x[c,p])        
    ss.x[c,p] ~ dunif(0,1) 
    ss.mu[c,p] <- logit(ss.x[c,p])
    f.x[c,p] ~ dunif(0,0.1) 
    f.mu[c,p] <- logit(f.x[c,p])               
    
    sa.sd[c,p] ~ dunif(0.05,2)                 # Priors for SDs of survival and recovery rates
    ss.sd[c,p] ~ dunif(0.05,2) 
    f.sd[c,p] ~ dunif(0.05,2) 
    sa.tau[c,p] <- pow(sa.sd[c,p],-2)          # Express as precision
    ss.tau[c,p] <- pow(ss.sd[c,p],-2) 
    f.tau[c,p] <- pow(f.sd[c,p],-2) 
    
    ## Generate annual parameter estimates
    for (t in 1:yrs){
    eps.sa[t,c,p] ~ dnorm(0,sa.tau[c,p])
    eps.ss[t,c,p] ~ dnorm(0,ss.tau[c,p])
    eps.f[t,c,p] ~ dnorm(0,f.tau[c,p])
    
    logit(sa[t,c,p]) <- sa.mu[c,p] + eps.sa[t,c,p]         # annual survival (by year, cohort, pop)
    logit(ss[t,c,p]) <- ss.mu[c,p] + eps.ss[t,c,p]         # seasonal survival (summer)
    logit(f[t,c,p]) <- f.mu[c,p] + eps.f[t,c,p]            # Brownie recovery rate
    
    } # close t
    
    ## Balance equations, assuming estimating N in both seasons
    # Initial pop sizes in summer [t, s, c, p]
    N[1,2,c,p] <- trunc(N[1,1,c,p]*ss[1,c,p])             # pop size in late summer is function of summer survival 
    } #c
    
    for (t in 2:yrs){
    N[t,1,1,p] <- trunc(N[t,1,3,p] * F[t,p])              # pop size of juv in spring = adF * mean fecundity
    N[t,2,1,p] <- trunc(N[t,1,1,p]*ss[t,1,p])             # pop size of juv in late summer
    sw[t,1,p] <- sa[t-1,1,p]/ss[t,1,p]                    # derived winter survival for juveniles
    
    for (c in 2:3){
    sw[t,c,p] <- sa[t-1,c,p]/ss[t,c,p]                                                      # derived winter survival for adults
    N[t,1,c,p] <- trunc(N[t-1,2,c,p]*sw[t,c,p] + pi.sex[t,c-1,p]*N[t-1,2,1,p]*sw[t,1,p])      # pop size of adults in spring
    N[t,2,c,p] <- trunc(N[t,1,c,p]*ss[t,c,p])                                               # pop size of adults in late summer 
    } #c
    } #t 
    
    # Observation process, recoveries and harvest data
    # Recoveries
    # Note: releases MUST be provided as data; we have saved as relAMWO
    for (c in 1:3){
    for (t in 1:yrs){
    marrayAMWO[t,1:(yrs+1),1,c,p] ~ dmulti(pr[t,,1,c,p], rel[t,1,c,p])        # Apr-Jun releases     
    marrayAMWO[t,1:(yrs+1),2,c,p] ~ dmulti(pr[t,,2,c,p], rel[t,2,c,p])        # Jul-Sep releases
    
    ## Define cell probabilities of m-arrays
    # Main diagonal--direct recovery in season [t]
    pr[t,t,1,c,p] <- ss[t,c,p] * f[t,c,p]                  # spring banded birds must survive to start of hunting season
    pr[t,t,2,c,p] <- f[t,c,p]                              # survival assumed 1 if banded Jul-Sep
    
    # monitor cumulative survival to start of next hunting season (previously surv; used in subsequent diagonals)
    cumS[t,t,1,c,p] <- ss[t,c,p] * sa[t,c,p]
    cumS[t,t,2,c,p] <- sa[t,c,p]
    } # t
    } # c
    
    for (t in 1:yrs){
    # Above main diagonal--indirect recovery in season [k > t]
    # All birds are adults now, no differences between Apr-Jun and Jul-Sep either
    for (k in (t+1):yrs){                                                                       # k loop to represent next year (above main diagonal)
    # recoveries
    pr[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi.sex[t,1,p]*f[k,2,p] + pi.sex[t,2,p]*f[k,3,p])          # juvs as mixture of AdM & AdF
    pr[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi.sex[t,1,p]*f[k,2,p] + pi.sex[t,2,p]*f[k,3,p])          # pi.sex[1] is proportion male; pi.sex[2] is proportion female      
    pr[t,k,1,2,p] <- cumS[t,k-1,1,2,p] * f[k,2,p] 
    pr[t,k,2,2,p] <- cumS[t,k-1,2,2,p] * f[k,2,p] 
    pr[t,k,1,3,p] <- cumS[t,k-1,1,3,p] * f[k,3,p] 
    pr[t,k,2,3,p] <- cumS[t,k-1,2,3,p] * f[k,3,p] 
    # monitor cumulative survival to start of next hunting period
    cumS[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi.sex[t,1,p]*sa[k,2,p] + pi.sex[t,2,p]*sa[k,3,p]) # juvs as mixture AdM & AdF
    cumS[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi.sex[t,1,p]*sa[k,2,p] + pi.sex[t,2,p]*sa[k,3,p]) 
    cumS[t,k,1,2,p] <- cumS[t,k-1,1,2,p] * sa[k,2,p]
    cumS[t,k,2,2,p] <- cumS[t,k-1,2,2,p] * sa[k,2,p]
    cumS[t,k,1,3,p] <- cumS[t,k-1,1,3,p] * sa[k,3,p]
    cumS[t,k,2,3,p] <- cumS[t,k-1,2,3,p] * sa[k,3,p]
    } #k
    
    # Left of main diag
    for (l in 1:(t-1)){           #l loop for previous year (left of main diagonal)                    
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
    
    # Bringing in harvest proportion data estimated from harvest model separately
    for (t in 1:yrs){ 
    for (c in 1:3){
    H[t,c,p] ~ dbin(h[t,c,p], N[t,2,c,p])
    
    h[t,c,p] <- f[t,c,p]/report[t]                    # harvest rate (h) is recovery rate divided by reporting rate (p)                                                      
    
    } #c                                              # note!! reporting rate now a vector of 0.506 and 0.8 when toll free bands started in 2001 (not just 0.506)
    } #t           
    } #p
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data
# Data to be pulled in from separate harvest model
H <- array(NA, dim = c(53,3,2))
pi <- array(NA, dim=c(53,3,2))

for(c in 1:3){
  for(p in 1:2){
    for(t in 1:51){
      pi[t,c,p] <- AMWO.harvest.wings.scaled.HPCC$mean$pi[t,c,p]            #pull in pi estimates
      H[t,c,p] <- round(AMWO.harvest.wings.scaled.HPCC$mean$H[t,c,p])       #pull in scaled harvest estimates (divided data by 10000 in previous model)
    }
  }
}

#replace missing H and pi values for last 2 years with means of last 5 years
H[52:53,1,] <- round(mean(H[47:51,1,], na.rm = TRUE))
H[52:53,2,] <- round(mean(H[47:51,2,], na.rm = TRUE))
H[52:53,3,] <- round(mean(H[47:51,3,], na.rm = TRUE))

pi[52:53,1,] <- mean(pi[47:51,1,], na.rm = TRUE)
pi[52:53,2,] <- mean(pi[47:51,2,], na.rm = TRUE)
pi[52:53,3,] <- mean(pi[47:51,3,], na.rm = TRUE)

#reporting rate vector. Toll free bands start in 2001 (using reporting rate of 0.8 for now, 0.506 for all years prior to that)
report<-c(0.506,rep(0.506,37),0.8,rep(0.8,14))

bugs.data <- list(yrs=dim(marrayAMWO)[1], marrayAMWO=marrayAMWO, rel=relAMWO, H=H, pi=pi, report=report)  

# inits
n.inits <- array(NA, dim=c(1,1,3,2))     #scale these inits by dividing by 10000
n.inits[1,1,3,1] <- 300
n.inits[1,1,3,2] <- 400
n.inits[1,1,2,1] <- 250
n.inits[1,1,2,2] <- 300

inits <- function(){list(
  n=n.inits
)}               

# Parameters monitored
parameters <- c("sa", "f", "N","sa.x", "ss.x", "f.x", "sa.sd", "ss.sd", "f.sd")  # anything else to include? h?

# MCMC settings
ni <- 350000
nt <- 10
nb <- 300000
nc <- 3

# Run JAGS
AMWO.combo.jags <- jagsUI(bugs.data, inits=inits, parameters, "Lincoln.Brownie.Simple2.6April.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)
save(AMWO.combo.jags, file="AMWO_Lincoln_Brownie_Combo_6May2017.Rda")
#load(file="AMWO_Lincoln_Brownie_Combo_6May2017.Rda")
#print(AMWO.combo.jags$summary, digits=3)

#y <- c(AMWO.combo.jags$mean$N[1:53,1,2,1])
#x <- seq(1,53,by=1)
#plot(x,y)