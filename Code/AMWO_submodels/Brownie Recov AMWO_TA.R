##################################################################################################
# Based off of Todd's ABDU Brownie model
# 2-season band recovery model for AMWO (summer and fall/winter)
# Female and Male, Brownie
# Juvs transition to adults after first winter
# Authors: T Arnold, S Saunders, M Farr, J Ribeiro, S Rossman, A Sussman, A Wright
##################################################################################################
# Indexing: 
#t for cohort year
#k for recovery year (if NE to cohort year)
#s for season [1=Apr-Jun, 2=Jul-Sep]
#c for cohort [1=juvenile, 2=adult male, 3=adult female]
#p for population [1=eastern, 2=central]
#------------------------------------------------------------------------------------------------

library(jagsUI)

sink("AMWO.Brownie2")
cat("
    model {
    
    ## Priors and constraints
    # format: x is mean on real scale, mu on logit scale
    # sa is survival.annual, ss is survival.summer, f is Brownie recovery rate
    # Todd's note: Brownie models measure true survival: we can call it s, not phi!
    
    for (p in 1:2){
    pi[p] <- 0.5                      # **proportion of banded juveniles that are male--set to 0.5
    for (c in 1:3){                   # following your format here, but my inclination is to recognize g =1:6 groups
    sa.x[c,p] ~ dunif(0,1) 
    sa.mu[c,p] <- logit(sa.x[c,p])    # neat trick to create uniform(0,1) prior on logit scale
    ss.x[c,p] ~ dunif(0,1) 
    ss.mu[c,p] <- logit(ss.x[c,p])
    f.x[c,p] ~ dunif(0,1) 
    f.mu[c,p] <- logit(f.x[c,p])
    
    sa.sd[c,p] ~ dunif(0.05,2)        # Priors for SDs of survival and recovery rates
    ss.sd[c,p] ~ dunif(0.05,2) 
    f.sd[c,p] ~ dunif(0.05,2) 
    sa.tau[c,p] <- pow(sa.sd[c,p],-2)    # Express as precision
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
    
    # Note: in jags (but not in WinBUGS), releases MUST be provided as data; we have saved as relAMWO
    marrayAMWO[t,1:(yrs+1),1,c,p] ~ dmulti(pr[t,,1,c,p], rel[t,1,c,p])        # SS: Apr-Jun releases     
    marrayAMWO[t,1:(yrs+1),2,c,p] ~ dmulti(pr[t,,2,c,p], rel[t,2,c,p])        # SS: Jul-Sep releases
    
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

AMWO.Brownie2.jags <- jagsUI(bugs.data, inits=inits, parameters, "AMWO.Brownie2", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE) #changed from inits=NULL
print(AMWO.Brownie2.jags,digits=2)
save(AMWO.Brownie2.jags, file="AMWO Brownie.rda")

#hist(AMWO.Brownie2.jags$sims.list$pi[,1])
#hist(AMWO.Brownie2.jags$sims.list$pi[,2])

mean(AMWO.Brownie2.jags$mean$sa[,1,1]) #mean juv S eastern: 0.262 [0.21-0.31]
mean(AMWO.Brownie2.jags$mean$sa[,2,1]) #mean adM S eastern: 0.346 [0.24-0.41]
mean(AMWO.Brownie2.jags$mean$sa[,3,1]) #mean adF S eastern: 0.476 [0.39-0.56]

mean(AMWO.Brownie2.jags$mean$sa[,1,2]) #mean juv S central: 0.248 [0.21-0.30]
mean(AMWO.Brownie2.jags$mean$sa[,2,2]) #mean adM S central: 0.444 [0.39-0.49]
mean(AMWO.Brownie2.jags$mean$sa[,3,2]) #mean adF S central: 0.479 [0.43-0.53]
