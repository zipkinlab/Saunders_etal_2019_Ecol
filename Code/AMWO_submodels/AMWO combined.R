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
  Central[i,1] <- harvest[i,1]  # year for C region
  Eastern[i,1] <- harvest[(i+49),1]  # year for E region
  for (j in 2:12){
    Central[i,j] <- harvest[(i),(j+12)] #harvest column to local recovs column from dataset
    Eastern[i,j] <- harvest[(i+49),(j+12)] 
  } # close j
  } # close i

head(Central)   #some NAs in data; 12 columns, 49 rows each
tail(Eastern)

# specify bugs model
sink("Lincoln.bug")
cat("
    model {
  # priors for initial population sizes (mean first 5 years Lincoln, SD 100,000)
    af.N[1,1] ~ dnorm(720000,1E-10) #assuming these come from previous Lincoln analysis; adult female central (1)
    af.N[2,1] ~ dnorm(600000,1E-10)  #adult female eastern (2)
    am.N[1,1] ~ dnorm(615000,1E-10)   #adult male central
    am.N[2,1] ~ dnorm(530000,1E-10)   #adult male eastern

  # Priors and constraints for population means and variances
    for (i in 1:2){   # 1 Central, 2 Eastern pop
      f.af.x[i] ~ dunif(0.005,0.2)          # uniform prior for mean Brownie recovery
      f.af.mu[i] <- logit(f.af.x[i])   #adult female mean recovery rate?
      f.am.x[i] ~ dunif(0.005,0.2)          
      f.am.mu[i] <- logit(f.af.x[i])       #adult male mean recov rate
      f.juv.x[i] ~ dunif(0.005,0.2)          
      f.juv.mu[i] <- logit(f.af.x[i])       #juvenile mean recov rate
      s.loc.x[i] ~ dunif (0,1)         # prior for S from local to HY stage
      s.loc.mu[i] <- logit(s.loc.x[i])   #prior for S from local to HY stage, mean

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
      F.x[i] ~ dunif(0,4) # logical bounds on fecundity (F) with 4 egg clutch
      afr.x[i] ~ dunif(-1,1) # uniform prior for mean growth rate, ad females 
      amr.x[i] ~ dunif(-1,1) # ad males
      F.sd[i] ~ dunif(0.01,1) 
      afr.sd[i] ~ dunif(0.01,1) # uniform prior for ann variation in mean growth rate
      amr.sd[i] ~ dunif(0.01,1) 
      F.tau[i] <- pow(F.sd[i],-2)
      afr.tau[i] <- pow(afr.sd[i],-2)
      amr.tau[i] <- pow(amr.sd[i],-2)

# Generate annual parameter estimates
    for (y in 1:yrs){                                   #yrs is 49 here
    # generate recovery rates
      logit.f.af[i,y] ~ dnorm(f.af.mu[i],f.af.tau[i])  #uses mean and variance from above, adult female
      logit.f.am[i,y] ~ dnorm(f.am.mu[i],f.am.tau[i])  #uses mean and variance from above, adult male
      logit.f.juv[i,y] ~ dnorm(f.juv.mu[i],f.juv.tau[i])   #uses mean and variance from above, juvenile
      logit.s.loc[i,y] ~ dnorm(s.loc.mu[i],s.loc.tau[i])   #uses mean and variance from above, S from local to HY stage
      logit(f.af[i,y]) <- logit.f.af[i,y]
      logit(f.am[i,y]) <- logit.f.am[i,y]
      logit(f.juv[i,y]) <- logit.f.juv[i,y]
      logit(s.loc[i,y]) <- logit.s.loc[i,y]
      f.loc[i,y] <- s.loc[i,y]*f.juv[i,y]     #true recovery rates for locals a product of survival and recovery of juveniles?

    # generate process components
      af.r[i,y] ~ dnorm(afr.x[i],afr.tau[i])   # need F (fecundity, not recov!) for each year, but last lambda not used, af.r is mean growth rate from above
      am.r[i,y] ~ dnorm(amr.x[i],amr.tau[i])   #am.r is mean growth rate for adult males pop segment
      F[i,y] ~ dnorm(F.x[i],F.tau[i])      # Fecundity cannot be <0 or >4
      af.lambda[i,y] <- exp(af.r[i,y])     # convert r to lambda by exponentiating, r = log(lambda), converting from stochastic population growth rate to finite growth rate
      am.lambda[i,y] <- exp(am.r[i,y])

      } # close y

  # State process (estimate Fecundity in year 1)
    juv.N[i,1] <- trunc(af.N[i,1]*F[i,1])   #trunc truncates towards 0. juv pop size function of adult females and fecundity
    for (y in 2:yrs){
      af.N[i,y] <- trunc(af.N[i,y-1] * af.lambda[i,y-1])   #adult female pop size function of previous yr times lambda
      am.N[i,y] <- trunc(am.N[i,y-1] * am.lambda[i,y-1])  #adult male pop size function of previous yr times lambda
      juv.N[i,y] <- trunc(af.N[i,y] * F[i,y])
    } # close y
  } # close i  

# observation process, recoveries and harvest data
   for (y in 1:yrs){
      C.af.R[y] ~ dbin(f.af[1,y],C.af.B[y])  #central adult female recovery data, bin distribution function of adult female recoveries and central region females banded
      C.am.R[y] ~ dbin(f.am[1,y],C.am.B[y])  #central adult male recovery data, function of adult male recoveries and central males banded
      C.juv.R[y] ~ dbin(f.juv[1,y],C.juv.B[y]) #central juvenile recovery data, function of juv recoveries and central imms banded
      C.loc.R[y] ~ dbin(f.loc[1,y],C.loc.B[y])  #central local recovery data, function of local recoveries and central locals banded
      E.af.R[y] ~ dbin(f.af[2,y],E.af.B[y])  #ditto for eastern
      E.am.R[y] ~ dbin(f.am[2,y],E.am.B[y])   #ditto for eastern
      E.juv.R[y] ~ dbin(f.juv[2,y],E.juv.B[y])  #ditto for eastern
      E.loc.R[y] ~ dbin(f.loc[2,y],E.loc.B[y])   #ditto for eastern
   }
   for (y in 2:yrs){
      C.af.H[y] ~ dbin(f.af[1,y], af.N[1,y])  #ad female harvest data dist. binom as adult female pop size (af.N) and recovery rate (f.af)
      C.am.H[y] ~ dbin(f.am[1,y], am.N[1,y])   #ad male harvest data dist. binom as adult male pop size (am.N) and recovery rate (f.am)
      C.juv.H[y] ~ dbin(f.juv[1,y], juv.N[1,y]) #juv harvest data (imm in excel sheet) dist. as juv pop size and recovery rate
      E.af.H[y] ~ dbin(f.af[2,y], af.N[2,y])  #ditto for eastern
      E.am.H[y] ~ dbin(f.am[2,y], am.N[2,y])   #ditto for eastern
      E.juv.H[y] ~ dbin(f.juv[2,y], juv.N[2,y])  #ditto for eastern
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


