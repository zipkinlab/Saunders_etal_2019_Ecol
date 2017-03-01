#-----------------------------------------------------------------------------------------------
# Indexing: 
#t for cohort year
#k for recovery year (if NE to cohort year)
#s for season [1=Apr-Jun, 2=Jul-Sep]
#c for cohort [1=juvenile, 2=adult male, 3=adult female]
#p for population [1=eastern, 2=central]
#------------------------------------------------------------------------------------------------

# Make sure AMWO marray and AMWO rel are loaded into workspace

# LOAD HARVEST DATA #
harvest <- read.csv("AMWO harvest.csv",header=TRUE)

# read data into proper matrix for jags: needs to be in t,c,p format, right?
#bring in year
clean<-matrix(NA,nrow=98,ncol=1)
clean<-data.frame(clean)
clean[,1]<- harvest$Season

#bring in population
clean[,2]<- harvest$Pop

#bring in age classes of harvest
clean[,3]<-harvest$AF
clean[,4]<-harvest$AM
clean[,5]<-harvest$AU
clean[,6]<-harvest$IF
clean[,7]<-harvest$IM
clean[,8]<-harvest$IU

colnames(clean)<-c("Year","Pop","AdF","AdM","AdU","JuvF","JuvM","JuvU")

#convert clean to t,c,p matrix called Harvest to put in Lincoln estimation below

Years<-unique(clean$Year)
Years<-sort(Years)          #SS sorted
NYear<-length(Years)
Region<-unique(clean$Pop)
Region<-sort(Region)       #SS sorted
NRegion<-length(Region)

Harvest<-array(NA,dim=c(NYear,NClass,NRegion),
           dimnames =list(Years,c("Juvenile","Adult_Male","Adult_Female"),
                          c("Eastern","Central")))
# to be continued...


#------------#
#-BUGS Model-#
#------------#

sink("Lincoln.jags")
cat("
    model {
  
    # Initial pop sizes??
    af.N[1,1] ~ dnorm(720000,1E-10) #assuming these come from previous Lincoln analysis; adult female central (1)
    af.N[2,1] ~ dnorm(600000,1E-10)  #adult female eastern (2)
    am.N[1,1] ~ dnorm(615000,1E-10)   #adult male central
    am.N[2,1] ~ dnorm(530000,1E-10)   #adult male eastern

    # Priors and constraints for population means and variances
    for (p in 1:2){                       # 1 Eastern, 2 Central
    pi[p] <- 0.5                          # **proportion of banded juveniles that are male--set to 0.5
    for (c in 1:3){                       # following your format here, but my inclination is to recognize g =1:6 groups
    sa.x[c,p] ~ dunif(0,1) 
    sa.mu[c,p] <- logit(sa.x[c,p])        # neat trick to create uniform(0,1) prior on logit scale
    ss.x[c,p] ~ dunif(0,1) 
    ss.mu[c,p] <- logit(ss.x[c,p])
    f.x[c,p] ~ dunif(0,1) 
    f.mu[c,p] <- logit(f.x[c,p])

    sa.sd[c,p] ~ dunif(0.05,2)            # Priors for SDs of survival and recovery rates
    ss.sd[c,p] ~ dunif(0.05,2) 
    f.sd[c,p] ~ dunif(0.05,2) 
    sa.tau[c,p] <- pow(sa.sd[c,p],-2)     # Express as precision
    ss.tau[c,p] <- pow(ss.sd[c,p],-2) 
    f.tau[c,p] <- pow(f.sd[c,p],-2) 

    # priors for annual growth rates and fecundities
    F.x[p] ~ dunif(0, 4)                       # logical bounds on fecundity (F) with 4 egg clutch; need c index??
    ar.x[2:3,p] ~ dunif(-1, 1)                 # uniform prior for mean growth rate, adults only? 
    F.sd[p] ~ dunif(0.01, 1) 
    ar.sd[2:3,p] ~ dunif(0.01, 1)              # uniform prior for ann variation in mean growth rate
    F.tau[p] <- pow(F.sd[p], -2)
    ar.tau[2:3,p] <- pow(ar.sd[2:3,p], -2)

    ## Generate annual parameter estimates
    for (t in 1:yrs){
    eps.sa[t,c,p] ~ dnorm(0,sa.tau[c,p])
    eps.ss[t,c,p] ~ dnorm(0,ss.tau[c,p])
    eps.f[t,c,p] ~ dnorm(0,f.tau[c,p])
    
    logit(sa[t,c,p]) <- sa.mu[c,p] + eps.sa[t,c,p]         # annual survival (by year, cohort, pop)
    logit(ss[t,c,p]) <- ss.mu[c,p] + eps.ss[t,c,p]         # seasonal survival (summer)
    logit(f[t,c,p]) <- f.mu[c,p] + eps.f[t,c,p]            # Brownie recovery rate

    # Generate process components
    a.r[t,2:3,p] ~ dnorm(ar.x[2:3,p], ar.tau[2:3,p])       # a.r is mean growth rate from above
    F[t,p] ~ dnorm(F.x[p], F.tau[p])                       # Fecundity cannot be <0 or >4
    a.lambda[t,2:3,p] <- exp(a.r[t,2:3,p])                 # convert r to lambda by exponentiating, r = log(lambda), converting from stochastic population growth rate to finite growth rate
    } # close t

    # State process (estimate Fecundity in year 1)
    N[1,1,p] <- trunc(N[1,3,p]*F[1,p])                          #trunc truncates towards 0. juv pop size function of adult females and fecundity
    for (t in 2:yrs){
      N[t,1,p] <- trunc(N[t,3,p] * F[t,p])
      N[t,2:3,p] <- trunc(N[t-1,2:3,p] * a.lambda[t-1,2:3,p])   #adult pop sizes function of previous yr times lambda
    } # close t

    # Observation process, recoveries and harvest data
    # Recoveries
    # Note: releases MUST be provided as data; we have saved as relAMWO
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
    }} # t,c
    
    for (t in 1:yrs){
    # Above main diagonal--indirect recovery in season [k > t]
    # All birds are adults now, no differences between Apr-Jun and Jul-Sep either
    for (k in (t+1):yrs){                                                                 # SS: k loop to represent next year (above main diagonal)
    # recoveries
    pr[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi[p]*f[k,2,p] + (1-pi[p])*f[k,3,p])            # juvs as mixture of AdM & AdF
    pr[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi[p]*f[k,2,p] + (1-pi[p])*f[k,3,p])          
    pr[t,k,1,2,p] <- cumS[t,k-1,1,2,p] * f[k,2,p] 
    pr[t,k,2,2,p] <- cumS[t,k-1,2,2,p] * f[k,2,p] 
    pr[t,k,1,3,p] <- cumS[t,k-1,1,3,p] * f[k,3,p] 
    pr[t,k,2,3,p] <- cumS[t,k-1,2,3,p] * f[k,3,p] 
    # monitor cumulative survival to start of next hunting period
    cumS[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi[p]*sa[k,2,p] + (1-pi[p])*sa[k,3,p])        # juvs as mixture AdM & AdF
    cumS[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi[p]*sa[k,2,p] + (1-pi[p])*sa[k,3,p]) 
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
    } #p

    for (t in 2:yrs){
    Harvest[t,c,p] ~ dmulti(f[t,c,p], N[t,c,p])           # harvest data function of pop size (N) and recovery rate (f)
    } #t
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(yrs=dim(marrayAMWO)[1], marrayAMWO=marrayAMWO, rel=relAMWO)

# Initial values (dimensions don't match 2 x 3 structure, so I'm using inits=NULL)       # SS: shouldn't be 3 x 2 structure? Made below, though
sa.x.inits <- matrix(NA, nrow=3, ncol=2)                                                 # may not be necessary.
fill <- runif(1,0,1)
for (p in 1:2){
  for (c in 1:3){
    sa.x.inits[c,p] <- fill
  }}

ss.x.inits <- matrix(NA, nrow=3, ncol=2)
fill <- runif(1,0,1)
for (p in 1:2){
  for (c in 1:3){
    ss.x.inits[c,p] <- fill
  }}

f.x.inits <- matrix(NA, nrow=3, ncol=2)
fill <- runif(1,0,1)
for (p in 1:2){
  for (c in 1:3){
    f.x.inits[c,p] <- fill
  }}

sa.sd.inits <- matrix(NA, nrow=3, ncol=2)
fill <- runif(1,0.05,2)
for (p in 1:2){
  for (c in 1:3){
    sa.sd.inits[c,p] <- fill
  }}

ss.sd.inits <- matrix(NA, nrow=3, ncol=2)
fill <- runif(1,0.05,2)
for (p in 1:2){
  for (c in 1:3){
    ss.sd.inits[c,p] <- fill
  }}

f.sd.inits <- matrix(NA, nrow=3, ncol=2)
fill <- runif(1,0.05,2)
for (p in 1:2){
  for (c in 1:3){
    f.sd.inits[c,p] <- fill
  }}

inits <- function(){list(sa.x=sa.x.inits, sa.sd=sa.sd.inits,
                         ss.x=ss.x.inits, ss.sd=ss.sd.inits,
                         f.x=f.x.inits, f.sd=f.sd.inits)}  

# Parameters monitored
parameters <- c("sa.x", "ss.x", "f.x", "sa.sd", "ss.sd", "f.sd", "sa", "f")        # no need to monitor pi if set to 0.5

# MCMC settings (everything converges)
ni <- 8000
nt <- 2
nb <- 2000
nc <- 3


