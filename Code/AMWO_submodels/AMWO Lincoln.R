#--------------------------------------------------------------------------------------------------------------
# This is some of Todd's earlier code that treats the Lincoln estimation as two stages
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
print(harvest[,1:2])
print(harvest[1,])

# partition data sets by region (East, Central) and time period (1964-1996, 1999-2013)
C.DSS <- matrix(0,nrow=33,ncol=12,byrow=TRUE)
C.HIP <- matrix(0,nrow=15,ncol=12,byrow=TRUE)
E.DSS <- matrix(0,nrow=33,ncol=12,byrow=TRUE)
E.HIP <- matrix(0,nrow=15,ncol=12,byrow=TRUE)

for (i in 1:33){
  C.DSS[i,1] <- harvest[i+1,1]  # year
  E.DSS[i,1] <- harvest[i+50,1]  
  for (j in 2:12){
    C.DSS[i,j] <- harvest[(i+1),(j+12)] 
    E.DSS[i,j] <- harvest[(i+50),(j+12)]} # close j
} # close i

for (i in 1:15){
  C.HIP[i,1] <- harvest[i+34,1]  # year
  E.HIP[i,1] <- harvest[i+83,1]  
  for (j in 2:12){
    C.HIP[i,j] <- harvest[i+34,j+12] 
    E.HIP[i,j] <- harvest[i+83,j+12]} # close j
} # close i

# inspect data
print(C.DSS)
head(E.DSS)
head(C.HIP)
head(E.HIP)

# specify bugs model
sink("Lincoln.bug")
cat("
    model {
    
    # Priors and constraints, for logit link on f
  for (i in 1:4){   # 1 Central DSS, 2 Eastern DSS, 3 Central HIP, 4 Eastern HIP
    f.af.x[i] ~ dunif(0,0.2)          # uniform prior for mean Brownie recovery, adult Female
    f.af.mu[i] <- logit(f.af.x[i])
    f.am.x[i] ~ dunif(0,0.2)          
    f.am.mu[i] <- logit(f.af.x[i])
    f.juv.x[i] ~ dunif(0,0.2)          
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
} # close group loop

# Generate annual parameter estimates for DSS data
  for (i in 1:2){
    for (y in 1:DSS){
    logit.f.af[i,y] ~ dnorm(f.af.mu[i],f.af.tau[i])
    logit.f.am[i,y] ~ dnorm(f.am.mu[i],f.am.tau[i])
    logit.f.juv[i,y] ~ dnorm(f.juv.mu[i],f.juv.tau[i])
    logit.s.loc[i,y] ~ dnorm(s.loc.mu[i],s.loc.tau[i])
    logit(f.af[i,y]) <- logit.f.af[i,y]
    logit(f.am[i,y]) <- logit.f.am[i,y]
    logit(f.juv[i,y]) <- logit.f.juv[i,y]
    logit(s.loc[i,y]) <- logit.s.loc[i,y]
    f.loc[i,y] <- s.loc[i,y]*f.juv[i,y]}
    }
# And for HIP
  for (i in 3:4){
    for (y in 1:HIP){
    logit.f.af[i,y] ~ dnorm(f.af.mu[i],f.af.tau[i])
    logit.f.am[i,y] ~ dnorm(f.am.mu[i],f.am.tau[i])
    logit.f.juv[i,y] ~ dnorm(f.juv.mu[i],f.juv.tau[i])
    logit.s.loc[i,y] ~ dnorm(s.loc.mu[i],s.loc.tau[i])
    logit(f.af[i,y]) <- logit.f.af[i,y]
    logit(f.am[i,y]) <- logit.f.am[i,y]
    logit(f.juv[i,y]) <- logit.f.juv[i,y]
    logit(s.loc[i,y]) <- logit.s.loc[i,y]
    f.loc[i,y] <- s.loc[i,y]*f.juv[i,y]}
    }
    
# Generate direct recovery rates and Lincoln estimates
  for (y in 1:DSS){
    C.DSS.af.R[y] ~ dbin(f.af[1,y],C.DSS.af.B[y])
    E.DSS.af.R[y] ~ dbin(f.af[2,y],E.DSS.af.B[y])
    C.DSS.am.R[y] ~ dbin(f.am[1,y],C.DSS.am.B[y])
    E.DSS.am.R[y] ~ dbin(f.am[2,y],E.DSS.am.B[y])
    C.DSS.juv.R[y] ~ dbin(f.juv[1,y],C.DSS.juv.B[y])
    E.DSS.juv.R[y] ~ dbin(f.juv[2,y],E.DSS.juv.B[y])
    C.DSS.loc.R[y] ~ dbin(f.loc[1,y],C.DSS.loc.B[y])
    E.DSS.loc.R[y] ~ dbin(f.loc[2,y],E.DSS.loc.B[y])

    C.DSS.af.N[y] <- C.DSS.af.H[y]/f.af[1,y]-1
    E.DSS.af.N[y] <- E.DSS.af.H[y]/f.af[2,y]-1
    C.DSS.am.N[y] <- C.DSS.am.H[y]/f.am[1,y]-1
    E.DSS.am.N[y] <- E.DSS.am.H[y]/f.am[2,y]-1
    C.DSS.juv.N[y] <- C.DSS.juv.H[y]/f.juv[1,y]-1
    E.DSS.juv.N[y] <- E.DSS.juv.H[y]/f.juv[2,y]-1
    }
  for (y in 1:HIP){
    C.HIP.af.R[y] ~ dbin(f.af[3,y],C.HIP.af.B[y])
    E.HIP.af.R[y] ~ dbin(f.af[4,y],E.HIP.af.B[y])
    C.HIP.am.R[y] ~ dbin(f.am[3,y],C.HIP.am.B[y])
    E.HIP.am.R[y] ~ dbin(f.am[4,y],E.HIP.am.B[y])
    C.HIP.juv.R[y] ~ dbin(f.juv[3,y],C.HIP.juv.B[y])
    E.HIP.juv.R[y] ~ dbin(f.juv[4,y],E.HIP.juv.B[y])
    C.HIP.loc.R[y] ~ dbin(f.loc[3,y],C.HIP.loc.B[y])
    E.HIP.loc.R[y] ~ dbin(f.loc[4,y],E.HIP.loc.B[y])

    C.HIP.af.N[y] <- C.HIP.af.H[y]/f.af[3,y]-1
    E.HIP.af.N[y] <- E.HIP.af.H[y]/f.af[4,y]-1
    C.HIP.am.N[y] <- C.HIP.am.H[y]/f.am[3,y]-1
    E.HIP.am.N[y] <- E.HIP.am.H[y]/f.am[4,y]-1
    C.HIP.juv.N[y] <- C.HIP.juv.H[y]/f.juv[3,y]-1
    E.HIP.juv.N[y] <- E.HIP.juv.H[y]/f.juv[4,y]-1
    }
    
  } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data  , add 1 to harvest, banding and recovery totals, per Seber's small-sample modification
bugs.data <- list(DSS= 33, HIP= 15, 
                  C.DSS.af.H=C.DSS[,2]+1, E.DSS.af.H=E.DSS[,2]+1, C.HIP.af.H=C.HIP[,2]+1, E.HIP.af.H=E.HIP[,2]+1,
                  C.DSS.am.H=C.DSS[,3]+1, E.DSS.am.H=E.DSS[,3]+1, C.HIP.am.H=C.HIP[,3]+1, E.HIP.am.H=E.HIP[,3]+1,
                  C.DSS.juv.H=C.DSS[,4]+1, E.DSS.juv.H=E.DSS[,4]+1, C.HIP.juv.H=C.HIP[,4]+1, E.HIP.juv.H=E.HIP[,4]+1,
                  C.DSS.af.B=C.DSS[,5]+1, E.DSS.af.B=E.DSS[,5]+1, C.HIP.af.B=C.HIP[,5]+1, E.HIP.af.B=E.HIP[,5]+1,
                  C.DSS.am.B=C.DSS[,6]+1, E.DSS.am.B=E.DSS[,6]+1, C.HIP.am.B=C.HIP[,6]+1, E.HIP.am.B=E.HIP[,6]+1,
                  C.DSS.juv.B=C.DSS[,7]+1, E.DSS.juv.B=E.DSS[,7]+1, C.HIP.juv.B=C.HIP[,7]+1, E.HIP.juv.B=E.HIP[,7]+1,
                  C.DSS.loc.B=C.DSS[,8]+1, E.DSS.loc.B=E.DSS[,8]+1, C.HIP.loc.B=C.HIP[,8]+1, E.HIP.loc.B=E.HIP[,8]+1,
                  C.DSS.af.R=C.DSS[,9]+1, E.DSS.af.R=E.DSS[,9]+1, C.HIP.af.R=C.HIP[,9]+1, E.HIP.af.R=E.HIP[,9]+1,
                  C.DSS.am.R=C.DSS[,10]+1, E.DSS.am.R=E.DSS[,10]+1, C.HIP.am.R=C.HIP[,10]+1, E.HIP.am.R=E.HIP[,10]+1,
                  C.DSS.juv.R=C.DSS[,11]+1, E.DSS.juv.R=E.DSS[,11]+1, C.HIP.juv.R=C.HIP[,11]+1, E.HIP.juv.R=E.HIP[,11]+1,
                  C.DSS.loc.R=C.DSS[,12]+1, E.DSS.loc.R=E.DSS[,12]+1, C.HIP.loc.R=C.HIP[,12]+1, E.HIP.loc.R=E.HIP[,12]+1)

# Initial values (not all priors listed)
inits <- function(){list(f.af.x=c(rep(runif(0,0.2),4)))}  

# Parameters monitored
parameters <- c("f.af.x", "f.am.x", "f.juv.x", "s.loc.x", 
                "C.DSS.af.N","E.DSS.af.N","C.HIP.af.N","E.HIP.af.N",
                "C.DSS.am.N","E.DSS.am.N","C.HIP.am.N","E.HIP.am.N",
                "C.DSS.juv.N","E.DSS.juv.N","C.HIP.juv.N","E.HIP.juv.N")

# MCMC settings 
ni <- 150000
nt <- 10        
nb <- 50000
nc <- 3

ABDU.Lincoln.jags <- jagsUI(bugs.data, inits, parameters, "Lincoln.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(ABDU.Lincoln.jags, digits = 1)
names(ABDU.Lincoln.jags$mean$C.DSS.af.N)

# create matrix to store Lincoln estimates and SD as output for part 2
Lincoln.DSS.out <- matrix(NA, nrow = 33, ncol = 13) 
Lincoln.HIP.out <- matrix(NA, nrow = 15, ncol = 13) 

# save Lincoln output to data file (truncate N, convert SD to tau)
for (y in 1:33){
  Lincoln.DSS.out[y,1]=C.DSS[y,1]
  Lincoln.DSS.out[y,2]=trunc(ABDU.Lincoln.jags$mean$C.DSS.af.N[y])
  Lincoln.DSS.out[y,3]=trunc(ABDU.Lincoln.jags$mean$E.DSS.af.N[y])
  Lincoln.DSS.out[y,4]=trunc(ABDU.Lincoln.jags$mean$C.DSS.am.N[y])
  Lincoln.DSS.out[y,5]=trunc(ABDU.Lincoln.jags$mean$E.DSS.am.N[y])
  Lincoln.DSS.out[y,6]=trunc(ABDU.Lincoln.jags$mean$C.DSS.juv.N[y])
  Lincoln.DSS.out[y,7]=trunc(ABDU.Lincoln.jags$mean$E.DSS.juv.N[y])
  Lincoln.DSS.out[y,8]=1/(ABDU.Lincoln.jags$sd$C.DSS.af.N[y])^2
  Lincoln.DSS.out[y,9]=1/(ABDU.Lincoln.jags$sd$E.DSS.af.N[y])^2
  Lincoln.DSS.out[y,10]=1/(ABDU.Lincoln.jags$sd$C.DSS.am.N[y])^2
  Lincoln.DSS.out[y,11]=1/(ABDU.Lincoln.jags$sd$E.DSS.am.N[y])^2
  Lincoln.DSS.out[y,12]=1/(ABDU.Lincoln.jags$sd$C.DSS.juv.N[y])^2
  Lincoln.DSS.out[y,13]=1/(ABDU.Lincoln.jags$sd$E.DSS.juv.N[y])^2
}

for (y in 1:15){
  Lincoln.HIP.out[y,1]=C.HIP[y,1]
  Lincoln.HIP.out[y,2]=trunc(ABDU.Lincoln.jags$mean$C.HIP.af.N[y])
  Lincoln.HIP.out[y,3]=trunc(ABDU.Lincoln.jags$mean$E.HIP.af.N[y])
  Lincoln.HIP.out[y,4]=trunc(ABDU.Lincoln.jags$mean$C.HIP.am.N[y])
  Lincoln.HIP.out[y,5]=trunc(ABDU.Lincoln.jags$mean$E.HIP.am.N[y])
  Lincoln.HIP.out[y,6]=trunc(ABDU.Lincoln.jags$mean$C.HIP.juv.N[y])
  Lincoln.HIP.out[y,7]=trunc(ABDU.Lincoln.jags$mean$E.HIP.juv.N[y])
  Lincoln.HIP.out[y,8]=1/(ABDU.Lincoln.jags$sd$C.HIP.af.N[y])^2
  Lincoln.HIP.out[y,9]=1/(ABDU.Lincoln.jags$sd$E.HIP.af.N[y])^2
  Lincoln.HIP.out[y,10]=1/(ABDU.Lincoln.jags$sd$C.HIP.am.N[y])^2
  Lincoln.HIP.out[y,11]=1/(ABDU.Lincoln.jags$sd$E.HIP.am.N[y])^2
  Lincoln.HIP.out[y,12]=1/(ABDU.Lincoln.jags$sd$C.HIP.juv.N[y])^2
  Lincoln.HIP.out[y,13]=1/(ABDU.Lincoln.jags$sd$E.HIP.juv.N[y])^2
}
head(Lincoln.DSS.out)
head(Lincoln.Hip.out)
tail(Lincoln.HIP.out)
write.csv(Lincoln.DSS.out,file="Lincoln DSS.csv")
write.csv(Lincoln.HIP.out,file="Lincoln HIP.csv")


#################################################
# Part 2: simple state-space IPM on Lincoln data
# log-normal version
################################################

sink("AMWO.SS.bug")  
cat("
    model {
    # priors for initial population sizes (too many, read data from file)
    C.DSS.af.N[1] ~ dnorm(C.DSS.af.mu,C.DSS.af.tau) # observed initial Lincoln estimates and precision (read in as data)
    E.DSS.af.N[1] ~ dnorm(E.DSS.af.mu,E.DSS.af.tau)
    C.DSS.am.N[1] ~ dnorm(C.DSS.am.mu,C.DSS.am.tau)
    E.DSS.am.N[1] ~ dnorm(E.DSS.am.mu,E.DSS.am.tau)

    C.HIP.af.N[1] ~ dnorm(C.HIP.af.mu,C.HIP.af.tau) # observed initial Lincoln estimates and precision (read in as data)
    E.HIP.af.N[1] ~ dnorm(E.HIP.af.mu,E.HIP.af.tau)
    C.HIP.am.N[1] ~ dnorm(C.HIP.am.mu,C.HIP.am.tau)
    E.HIP.am.N[1] ~ dnorm(E.HIP.am.mu,E.HIP.am.tau)

  # priors for annual growth rates and fecundities
    for (i in 1:4){
      F.x[i] ~ dunif(0,4) # C.DSS, E.DSS, C.HIP, E.HIP
      r.x[i] ~ dunif(-1,1) # uniform prior for mean growth rate C.DSS.af, E.DSS.af, C.DSS.am, E.DSS.am
      r.x[i+4] ~ dunif(-1,1) # HIP series
      F.sd[i] ~ dunif(0.01,1) # C.DSS, E.DSS, C.HIP, E.HIP
      r.sd[i] ~ dunif(0.01,1) # uniform prior for ann variation in mean growth rate
      r.sd[i+4] ~ dunif(0.01,1) # HIP series
      F.tau[i] <- pow(F.sd[i],-2)
      r.tau[i] <- pow(r.sd[i],-2)
      r.tau[i+4] <- pow(r.sd[i+4],-2)
    } # close i loop
      for (y in 1:DSS){
        C.DSS.af.r[y] ~ dnorm(r.x[1],r.tau[1])   # need F for each year, but last lambda will not be used
        E.DSS.af.r[y] ~ dnorm(r.x[2],r.tau[2])
        C.DSS.am.r[y] ~ dnorm(r.x[3],r.tau[3])
        E.DSS.am.r[y] ~ dnorm(r.x[4],r.tau[4])
        C.DSS.F[y] ~ dnorm(F.x[1],F.tau[1])      # Fecundity cannot be <0 or >4
        E.DSS.F[y] ~ dnorm(F.x[2],F.tau[2])
        C.DSS.af.lambda[y] <- exp(C.DSS.af.r[y]) # convert r to lambda by exponentiating, r = log(lambda)
        E.DSS.af.lambda[y] <- exp(E.DSS.af.r[y])
        C.DSS.am.lambda[y] <- exp(C.DSS.am.r[y])
        E.DSS.am.lambda[y] <- exp(E.DSS.am.r[y])
      } # close y loop for DSS
      for (y in 1:HIP){
        C.HIP.af.r[y] ~ dnorm(r.x[5],r.tau[5])
        E.HIP.af.r[y] ~ dnorm(r.x[6],r.tau[6])
        C.HIP.am.r[y] ~ dnorm(r.x[7],r.tau[7])
        E.HIP.am.r[y] ~ dnorm(r.x[8],r.tau[8])
        C.HIP.F[y] ~ dnorm(F.x[3],F.tau[3])
        E.HIP.F[y] ~ dnorm(F.x[4],F.tau[4])
        C.HIP.af.lambda[y] <- exp(C.HIP.af.r[y])
        E.HIP.af.lambda[y] <- exp(E.HIP.af.r[y])
        C.HIP.am.lambda[y] <- exp(C.HIP.am.r[y])
        E.HIP.am.lambda[y] <- exp(E.HIP.am.r[y])
      } # close y loop for HIP

    # Likelihood
    # State process (estimate Fecundity in year 1)
    C.DSS.juv.N[1] <- C.DSS.af.N[1]*C.DSS.F[1]
    E.DSS.juv.N[1] <- E.DSS.af.N[1]*E.DSS.F[1]
    C.HIP.juv.N[1] <- C.HIP.af.N[1]*C.HIP.F[1]
    E.HIP.juv.N[1] <- C.HIP.af.N[1]*E.HIP.F[1]

    for (y in 2:DSS){
    C.DSS.af.N[y] <- C.DSS.af.N[y-1] * C.DSS.af.lambda[y-1]
    E.DSS.af.N[y] <- E.DSS.af.N[y-1] * E.DSS.af.lambda[y-1]
    C.DSS.am.N[y] <- C.DSS.am.N[y-1] * C.DSS.am.lambda[y-1]
    E.DSS.am.N[y] <- E.DSS.am.N[y-1] * E.DSS.am.lambda[y-1]
    C.DSS.juv.N[y] <- C.DSS.af.N[y] * C.DSS.F[y]
    E.DSS.juv.N[y] <- E.DSS.af.N[y] * E.DSS.F[y]} # close y
      

    for (y in 2:HIP){
    C.HIP.af.N[y] <- C.HIP.af.N[y-1] * C.HIP.af.lambda[y-1]
    E.HIP.af.N[y] <- E.HIP.af.N[y-1] * E.HIP.af.lambda[y-1]
    C.HIP.am.N[y] <- C.HIP.am.N[y-1] * C.HIP.am.lambda[y-1]
    E.HIP.am.N[y] <- E.HIP.am.N[y-1] * E.HIP.am.lambda[y-1]
    C.HIP.juv.N[y] <- C.HIP.af.N[y] * C.HIP.F[y]
    E.HIP.juv.N[y] <- E.HIP.af.N[y] * E.HIP.F[y]} # close y

    # Observation process
    for (y in 1:DSS) {
    C.DSS.af.L[y] ~ dnorm(C.DSS.af.N[y],C.DSS.af.prec[y])
    E.DSS.af.L[y] ~ dnorm(E.DSS.af.N[y],E.DSS.af.prec[y])
    C.DSS.am.L[y] ~ dnorm(C.DSS.am.N[y],C.DSS.am.prec[y])
    E.DSS.am.L[y] ~ dnorm(E.DSS.am.N[y],E.DSS.am.prec[y])
    C.DSS.juv.L[y] ~ dnorm(C.DSS.juv.N[y],C.DSS.juv.prec[y])
    E.DSS.juv.L[y] ~ dnorm(E.DSS.juv.N[y],E.DSS.juv.prec[y])
    }
    for (y in 1:HIP) {
    C.HIP.af.L[y] ~ dnorm(C.HIP.af.N[y],C.HIP.af.prec[y])
    E.HIP.af.L[y] ~ dnorm(E.HIP.af.N[y],E.HIP.af.prec[y])
    C.HIP.am.L[y] ~ dnorm(C.HIP.am.N[y],C.HIP.am.prec[y])
    E.HIP.am.L[y] ~ dnorm(E.HIP.am.N[y],E.HIP.am.prec[y])
    C.HIP.juv.L[y] ~ dnorm(C.HIP.juv.N[y],C.HIP.juv.prec[y])
    E.HIP.juv.L[y] ~ dnorm(E.HIP.juv.N[y],E.HIP.juv.prec[y])
    }
    
    }    
    ",fill=TRUE)
sink()

# Bundle data  (swithched inits to median from initials because initials were outliers)
bugs.data <- list(DSS= 33, HIP= 15, 
                  C.DSS.af.mu=median(Lincoln.DSS.out[,2]), C.DSS.af.tau=median(Lincoln.DSS.out[,8]), 
                  E.DSS.af.mu=median(Lincoln.DSS.out[,3]), E.DSS.af.tau=median(Lincoln.DSS.out[,9]), 
                  C.DSS.am.mu=median(Lincoln.DSS.out[,4]), C.DSS.am.tau=median(Lincoln.DSS.out[,10]), 
                  E.DSS.am.mu=median(Lincoln.DSS.out[,5]), E.DSS.am.tau=median(Lincoln.DSS.out[,11]), 
                  C.HIP.af.mu=median(Lincoln.HIP.out[,2]), C.HIP.af.tau=median(Lincoln.HIP.out[,8]), 
                  E.HIP.af.mu=median(Lincoln.HIP.out[,3]), E.HIP.af.tau=median(Lincoln.HIP.out[,9]), 
                  C.HIP.am.mu=median(Lincoln.HIP.out[,4]), C.HIP.am.tau=median(Lincoln.HIP.out[,10]), 
                  E.HIP.am.mu=median(Lincoln.HIP.out[,5]), E.HIP.am.tau=median(Lincoln.HIP.out[,11]), 
                  C.DSS.af.L=Lincoln.DSS.out[,2], C.DSS.af.prec=Lincoln.DSS.out[,8], 
                  E.DSS.af.L=Lincoln.DSS.out[,3], E.DSS.af.prec=Lincoln.DSS.out[,9], 
                  C.DSS.am.L=Lincoln.DSS.out[,4], C.DSS.am.prec=Lincoln.DSS.out[,10], 
                  E.DSS.am.L=Lincoln.DSS.out[,5], E.DSS.am.prec=Lincoln.DSS.out[,11], 
                  C.DSS.juv.L=Lincoln.DSS.out[,6], C.DSS.juv.prec=Lincoln.DSS.out[,12], 
                  E.DSS.juv.L=Lincoln.DSS.out[,7], E.DSS.juv.prec=Lincoln.DSS.out[,13], 
                  C.HIP.af.L=Lincoln.HIP.out[,8], C.HIP.af.prec=Lincoln.HIP.out[,8], 
                  E.HIP.af.L=Lincoln.HIP.out[,9], E.HIP.af.prec=Lincoln.HIP.out[,9], 
                  C.HIP.am.L=Lincoln.HIP.out[,10], C.HIP.am.prec=Lincoln.HIP.out[,10], 
                  E.HIP.am.L=Lincoln.HIP.out[,11], E.HIP.am.prec=Lincoln.HIP.out[,11], 
                  C.HIP.juv.L=Lincoln.HIP.out[,12], C.HIP.juv.prec=Lincoln.HIP.out[,12], 
                  E.HIP.juv.L=Lincoln.HIP.out[,13], E.HIP.juv.prec=Lincoln.HIP.out[,13])

# Initial values  
inits <- function() {list(F.x=c(rep(runif(0,4),4)))}

# Parameters monitored   
parameters <- c("F.x", "r.x", "C.DSS.af.N", "E.DSS.af.N", "C.DSS.am.N", "E.DSS.am.N", "C.DSS.juv.N", "E.DSS.juv.N",
                "C.HIP.af.N", "E.HIP.af.N", "C.HIP.am.N", "E.HIP.am.N", "C.HIP.juv.N", "E.HIP.juv.N")

# MCMC settings (150000 w 50000 burn gives good convergence on all parms R-hat <= 1.01)
ni <- 120000
nt <- 10
nb <- 20000
nc <- 3

# jags version works well
AMWO.SS.jags <- jagsUI(bugs.data, inits, parameters, "AMWO.SS.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
print(AMWO.SS.jags,digits=3)
