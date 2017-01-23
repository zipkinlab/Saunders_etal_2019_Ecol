#--------------------------------------------------------------------------------------------------------------
# This is the integration of Lincoln estimation into the SS model itself via binomial distribution (i.e. Matt's way)
# Needs verification
# Based on simple dataset (AMWO harvest.csv)--can't decompose vital rates with this one
# Code sent 1/20/17
#--------------------------------------------------------------------------------------------------------------

# Lincoln Estimates for American Woodcock
# Todd Arnold, Univ. of Minnesota, arnol065@umn.edu
rm(list=ls(all=TRUE)) # clear memory

# load R2WinBUGS, set working directory for data and output, and identify WinBUGS location for your computer
library(jagsUI)
setwd("C:/Users/arnol065/Documents/Publications/AMWO IPM") 

# LOAD DATA #
harvest <- read.csv("AMWO harvest.csv",header=TRUE)

# read data into proper matrix for jags
Central <- matrix(0,nrow=49,ncol=12,byrow=TRUE)
Eastern <- matrix(0,nrow=49,ncol=12,byrow=TRUE)

for (i in 1:49){
  Central[i,1] <- harvest[i,1]  # year
  Eastern[i,1] <- harvest[(i+49),1]  # year
  for (j in 2:12){
    Central[i,j] <- harvest[(i),(j+12)] 
    Eastern[i,j] <- harvest[(i+49),(j+12)] 
  } # close j
  } # close i

head(Central)
tail(Eastern)

# specify bugs model
sink("Lincoln.bug")
cat("
    model {
  # priors for initial population sizes (mean first 5 years Lincoln, SD 100,000)
    af.N[1,1] ~ dnorm(720000,1E-10) 
    af.N[2,1] ~ dnorm(600000,1E-10)
    am.N[1,1] ~ dnorm(615000,1E-10) 
    am.N[2,1] ~ dnorm(530000,1E-10)

  # Priors and constraints for population means and variances
    for (i in 1:2){   # 1 Central, 2 Eastern pop
      f.af.x[i] ~ dunif(0.005,0.2)          # uniform prior for mean Brownie recovery
      f.af.mu[i] <- logit(f.af.x[i])
      f.am.x[i] ~ dunif(0.005,0.2)          
      f.am.mu[i] <- logit(f.af.x[i])
      f.juv.x[i] ~ dunif(0.005,0.2)          
      f.juv.mu[i] <- logit(f.af.x[i])
      s.loc.x[i] ~ dunif (0,1)         # prior for S from local to HY stage
      s.loc.mu[i] <- logit(s.loc.x[i])

    # priors for annual SD of annual recovery rate (logit scale)
      f.af.sd[i] ~ dunif(0.05,2)
      f.am.sd[i] ~ dunif(0.05,2)
      f.juv.sd[i] ~ dunif(0.05,2)
      s.loc.sd[i] ~ dunif (0.05,2)
      f.af.tau[i] <- pow(f.af.sd[i],-2)
      f.am.tau[i] <- pow(f.am.sd[i],-2)
      f.juv.tau[i] <- pow(f.juv.sd[i],-2)
      s.loc.tau[i] <- pow(s.loc.sd[i],-2)

  # priors for annual growth rates and fecundities
      F.x[i] ~ dunif(0,4) # logical bounds on fecundity with 4 egg clutch
      afr.x[i] ~ dunif(-1,1) # uniform prior for mean growth rate, ad females 
      amr.x[i] ~ dunif(-1,1) # ad males
      F.sd[i] ~ dunif(0.01,1) 
      afr.sd[i] ~ dunif(0.01,1) # uniform prior for ann variation in mean growth rate
      amr.sd[i] ~ dunif(0.01,1) 
      F.tau[i] <- pow(F.sd[i],-2)
      afr.tau[i] <- pow(afr.sd[i],-2)
      amr.tau[i] <- pow(amr.sd[i],-2)

# Generate annual parameter estimates
    for (y in 1:yrs){
    # generate recovery rates
      logit.f.af[i,y] ~ dnorm(f.af.mu[i],f.af.tau[i])
      logit.f.am[i,y] ~ dnorm(f.am.mu[i],f.am.tau[i])
      logit.f.juv[i,y] ~ dnorm(f.juv.mu[i],f.juv.tau[i])
      logit.s.loc[i,y] ~ dnorm(s.loc.mu[i],s.loc.tau[i])
      logit(f.af[i,y]) <- logit.f.af[i,y]
      logit(f.am[i,y]) <- logit.f.am[i,y]
      logit(f.juv[i,y]) <- logit.f.juv[i,y]
      logit(s.loc[i,y]) <- logit.s.loc[i,y]
      f.loc[i,y] <- s.loc[i,y]*f.juv[i,y]

    # generate process components
      af.r[i,y] ~ dnorm(afr.x[i],afr.tau[i])   # need F for each year, but last lambda not used
      am.r[i,y] ~ dnorm(amr.x[i],amr.tau[i])
      F[i,y] ~ dnorm(F.x[i],F.tau[i])      # Fecundity cannot be <0 or >4
      af.lambda[i,y] <- exp(af.r[i,y])     # convert r to lambda by exponentiating, r = log(lambda)
      am.lambda[i,y] <- exp(am.r[i,y])

      } # close y

  # State process (estimate Fecundity in year 1)
    juv.N[i,1] <- trunc(af.N[i,1]*F[i,1])
    for (y in 2:yrs){
      af.N[i,y] <- trunc(af.N[i,y-1] * af.lambda[i,y-1])
      am.N[i,y] <- trunc(am.N[i,y-1] * am.lambda[i,y-1])
      juv.N[i,y] <- trunc(af.N[i,y] * F[i,y])
    } # close y
  } # close i  

# observation process, recoveries and harvest data
   for (y in 1:yrs){
      C.af.R[y] ~ dbin(f.af[1,y],C.af.B[y])
      C.am.R[y] ~ dbin(f.am[1,y],C.am.B[y])
      C.juv.R[y] ~ dbin(f.juv[1,y],C.juv.B[y])
      C.loc.R[y] ~ dbin(f.loc[1,y],C.loc.B[y])
      E.af.R[y] ~ dbin(f.af[2,y],E.af.B[y])
      E.am.R[y] ~ dbin(f.am[2,y],E.am.B[y])
      E.juv.R[y] ~ dbin(f.juv[2,y],E.juv.B[y])
      E.loc.R[y] ~ dbin(f.loc[2,y],E.loc.B[y])
   }
   for (y in 2:yrs){
      C.af.H[y] ~ dbin(f.af[1,y], af.N[1,y])
      C.am.H[y] ~ dbin(f.am[1,y], am.N[1,y])
      C.juv.H[y] ~ dbin(f.juv[1,y], juv.N[1,y])
      E.af.H[y] ~ dbin(f.af[2,y], af.N[2,y])
      E.am.H[y] ~ dbin(f.am[2,y], am.N[2,y])
      E.juv.H[y] ~ dbin(f.juv[2,y], juv.N[2,y])
   }


  } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data  , add 1 to harvest, banding and recovery totals, per Seber's small-sample modification
bugs.data <- list(yrs = 49, C.af.H=trunc(Central[,2]+1), E.af.H=trunc(Eastern[,2]+1), 
                  C.am.H=trunc(Central[,3]+1), E.am.H=trunc(Eastern[,3]+1), 
                  C.juv.H=trunc(Central[,4]+1), E.juv.H=trunc(Eastern[,4]+1), 
                  C.af.B=Central[,5]+1, E.af.B=Eastern[,5]+1, 
                  C.am.B=Central[,6]+1, E.am.B=Eastern[,6]+1, 
                  C.juv.B=Central[,7]+1, E.juv.B=Eastern[,7]+1, 
                  C.loc.B=Central[,8]+1, E.loc.B=Eastern[,8]+1, 
                  C.af.R=Central[,9]+1, E.af.R=Eastern[,9]+1, 
                  C.am.R=Central[,10]+1, E.am.R=Eastern[,10]+1, 
                  C.juv.R=Central[,11]+1, E.juv.R=Eastern[,11]+1, 
                  C.loc.R=Central[,12]+1, E.loc.R=Eastern[,12]+1)
bugs.data$C.af.H
# Initial values (not all priors listed)
inits <- function(){list(f.af.x=c(rep(0.05,2)),f.am.x=c(rep(0.05,2)),f.juv.x=c(rep(0.05,2)))}  

# Parameters monitored
parameters <- c("F.x", "afr.x", "amr.x", "f.af.x", "f.am.x", "f.juv.x", "s.loc.x", 
                "af.N","am.N","juv.N")

# MCMC settings (just trying to get this to move)
ni <- 150000
nt <- 10        
nb <- 50000
nc <- 3

# hangs on first estimate
ABDU.Lincoln.jags <- jagsUI(bugs.data, inits, parameters, "Lincoln.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)

print(ABDU.Lincoln.jags, digits = 1)


