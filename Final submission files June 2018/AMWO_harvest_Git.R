#--------------------------------------------------------------------------------------------------------------
# American woodcock harvest model, part of AMWO IPM
# Code written by Saunders et al. January 2017 - December 2017
# Estimating annual harvest of AMWO during duck stamp surveys, harvest information program, and parts-collection survey
#--------------------------------------------------------------------------------------------------------------

library(jagsUI)

# Load data objects for HPCC run
load(file="Harvest_model_data.Rda")

#############################################################
# specify jags model to predict total annual harvest
#############################################################
sink("AMWO.state.harvest.pi.jags")
cat("
    model {
    
    # priors for state specific 2013 harvest (max of each component is 5-year mean total harvest)
    for (i in 1:n.states){
    Harvest.a[i,1] ~ dunif(0,mean5.HIP[i]) 
    Harvest.b[i,1] ~ dunif(0,mean5.HIP[i]) 
    Harvest.c[i,1] ~ dunif(0,mean5.HIP[i]) 
    Harvest[i,1] <- Harvest.a[i,1] + Harvest.b[i,1] + Harvest.c[i,1] # total harvest, derived parm
    
    # priors for latent process variation in components of true harvest
    sigma.Ha[i] ~ dunif(1,20000)  
    tau.Ha[i] <- pow(sigma.Ha[i],-2)
    sigma.Hb[i] ~ dunif(1,20000)
    tau.Hb[i] <- pow(sigma.Hb[i],-2)
    sigma.Hc[i] ~ dunif(1,20000)
    tau.Hc[i] <- pow(sigma.Hc[i],-2)
    
    }
    
    # State process (estimate true Harvest in years 2-51)
    for (i in 1:n.states){
    for (j in 2:n.yrs){
    eps.a[i,j] ~ dnorm(0,tau.Ha[i]) 
    eps.b[i,j] ~ dnorm(0,tau.Hb[i])
    eps.c[i,j] ~ dnorm(0,tau.Hc[i])
    Harvest.a[i,j] <- max(0,Harvest.a[i,j-1] + eps.a[i,j]) # autoregressive model, don't allow negative harvest
    Harvest.b[i,j] <- max(0,Harvest.b[i,j-1] + eps.b[i,j]) 
    Harvest.c[i,j] <- max(0,Harvest.c[i,j-1] + eps.c[i,j]) 
    Harvest[i,j] <- Harvest.a[i,j] + Harvest.b[i,j] + Harvest.c[i,j] # total harvest, derived parm
    }
    }
    
    # observation process, recoveries and harvest data
    for (i in 1:n.states){
    for (j in 1:n.HIP){
    HIP[i,j] ~ dnorm(Harvest[i,j],tau.HIP[i,j])} # end HIP years
    for (j in 13:n.yrs){
    frac[i,j] <- (stamps[i,j] - always.hunt[i])/(max.stamps[i] - always.hunt[i])
    DSS[i,j] ~ dnorm(Harvest.b[i,j]+frac[i,j]*Harvest.c[i,j],tau.DSS[i,j])} # end DSS yrs
    } # end i loop
    
    # calculating derived parameters of E and C harvest
    for (j in 1:n.yrs){
    Eastern.H[j] <- Harvest[2,j] + Harvest[3,j] + Harvest[4,j] + Harvest[5,j] + Harvest[6,j] + Harvest[7,j] + Harvest[11,j] + Harvest[13,j]
    Central.H[j] <- Harvest[1,j] + Harvest[8,j] + Harvest[9,j] + Harvest[10,j] + Harvest[12,j] + Harvest[14,j]
    
    region.H[j,1] <- Eastern.H[j]
    region.H[j,2] <- Central.H[j]
    }
    
    # Determining age-sex class proportions (pi's)
    for (p in 1:2){
    for (j in 1:n.yrs){
    pi[j,1,p] ~ dunif(0.1, 0.8)
    pi[j,2,p] ~ dunif(0.1, 0.8)
    pi[j,3,p] <- 1-(pi[j,1,p]+pi[j,2,p])
    
    # Summarize known wing samples by age and sex
    # populations (1 eastern, 2 central)
    wings.class[j,,p] ~ dmulti(pi[j,,p], wings.total[j,p])          # pi is 3 proportions: 1=juvs, 2=males, 3=females    
    H[j,1,p] <- pi[j,1,p]*region.H[rev[j],p]
    H[j,2,p] <- pi[j,2,p]*region.H[rev[j],p]
    H[j,3,p] <- pi[j,3,p]*region.H[rev[j],p]
    } #j
    } #p
    
    } # end jags model
    ",fill = TRUE)
sink()


# Bundle data, add 1 to harvest, banding and recovery totals, per Seber's small-sample modification
bugs.data <- list(n.yrs=n.yrs, n.states=n.states, n.HIP=n.HIP, mean5.HIP=mean5.HIP,
                  HIP=HIP, tau.HIP=tau.HIP, DSS=DSS, tau.DSS=tau.DSS,
                  stamps=stamps, max.stamps=max.stamps, always.hunt=always.hunt, rev=rev, wings.total=wings.total, wings.class=wings.class)

# Initial values (not all priors listed)
pi.inits <- array(dim=c(n.yrs,3,2))
for (c in 1:3){
  for (p in 1:2){
    for (j in 1:n.yrs){
      pi.inits[j,1,p] <- 0.46
      pi.inits[j,2,p] <- 0.21
      pi.inits[j,3,p] <- NA
    }}}

inits <- function(){list(pi=pi.inits)}

# Parameters monitored
parameters <- c("Harvest","Harvest.a","Harvest.b","Harvest.c", "Eastern.H","Central.H","pi","H")

# MCMC settings
na <- 1000
ni <- 150000
nt <- 10        
nb <- 50000
nc <- 3

# call jags
AMWO.state.harvest.pi <- jagsUI(bugs.data, inits, parameters, "AMWO.state.harvest.pi.jags", n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
save(AMWO.state.harvest.pi, file="AMWO_harvest_pi_Oct2017.Rda")
