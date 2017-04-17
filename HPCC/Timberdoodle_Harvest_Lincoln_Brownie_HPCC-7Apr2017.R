#--------------------------------------------------------------------------------------------------------------
# Estimating annual harvest of AMWO during duck stamp years using annual sales as a covariate
#--------------------------------------------------------------------------------------------------------------

#rm(list=ls(all=TRUE)) # clear memory

#load JAGS
library(jagsUI)

#load Data
load("Timberdoodle_HPCC_7Apr.R")

#########
## Part - MODEL 1
#########

bugs.data.1 <- list(yrs=51, HIP=HIP, DSS=DSS, stamps=stamps, wings.total=wings.total, wings.class=wings.class,
                  rev=rev)


# Initial values (not all priors listed)
pi.inits <- matrix(nrow=3, ncol=2)
for (c in 1:3){
  for (p in 1:2){
    pi.inits[1,p] <- 0.46
    pi.inits[2,p] <- 0.21
    pi.inits[3,p] <- NA
  }}

inits.1 <- function(){list(max.hunters=c(2200000,1000000), pi=pi.inits, harvest=runif(2,300000,500000))}  

# Parameters monitored
parameters.1 <- c("Harvest","max.hunters","sd.DSS","sd.HIP","sigma.H","H","pi")

# MCMC settings
ni.1 <- 1000
nt.1 <- 1        
#nb.1 <- 1500000
nc.1 <- 3

# call jags
AMWO.harvest.jags <- jagsUI(bugs.data.1, 
                            inits.1, 
                            parameters.1, 
                            "Timberdoodle_AMWO_harvest_wings.txt", 
                            n.chains = nc.1, 
                            n.thin = nt.1, 
                            n.iter = ni.1, 
                            #n.burnin = nb.1, 
                            parallel=TRUE, 
                            store.data=TRUE)

#save(AMWO.harvest.jags, file='AMWO_harvest_jags.R')

#########
## Part - MODEL 1
#########

# Bundle data

H <- array(NA, dim = c(53,3,2))
pi <- matrix(NA, nrow = 3, ncol = 2)

for(c in 1:3){
  for(p in 1:2){
    pi[c,p] <- AMWO.harvest.jags$mean$pi[c,p]
    for(t in 1:51){
      H[t,c,p] <- round(AMWO.harvest.jags$mean$H[t,c,p])
    }
  }
}

H[52:53,1,] <- round(mean(H[47:51,1,], na.rm = TRUE))
H[52:53,2,] <- round(mean(H[47:51,2,], na.rm = TRUE))
H[52:53,3,] <- round(mean(H[47:51,3,], na.rm = TRUE))

bugs.data.2 <- list(yrs=dim(marrayAMWO)[1], marrayAMWO=marrayAMWO, rel=relAMWO, H=H, pi=pi
                  )  

# inits
n.inits <- array(NA, dim=c(1,1,3,2))
n.inits[1,1,3,1] <- 2000000
n.inits[1,1,3,2] <- 1600000
n.inits[1,1,2,1] <- 1200000
n.inits[1,1,2,2] <- 1700000

inits.2 <- function(){list(
  n=n.inits
)}                 

# Parameters monitored
parameters.2 <- c("sa.x", "ss.x", "f.x", "sa.sd", "ss.sd", "f.sd",'F.x','F.sd', "sa", "f", "N")

# MCMC settings
ni.2 <- 1000
nt.2 <- 1        
#nb.1 <- 1500000
nc.2 <- 3

# Run JAGS
AMWO.combo.jags <- jagsUI(bugs.data.2, 
                          inits=inits.2, 
                          parameters.2, 
                          "Timberdoodle_Lincoln_Brownie_Simple2_6April.txt", 
                          n.chains = nc.2, 
                          n.thin = nt.2, 
                          n.iter = ni.2, 
                          #n.burnin = nb.2,
                          parallel=TRUE, 
                          store.data=TRUE) 

#save(AMWO.combo.jags, file="AMWO_combo_jags.R")




