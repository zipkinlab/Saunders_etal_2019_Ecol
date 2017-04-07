#--------------------------------------------------------------------------------------------------------------
# Estimating annual harvest of AMWO during duck stamp years using annual sales as a covariate
#--------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE)) # clear memory

# load R2WinBUGS, set working directory for data and output, and identify WinBUGS location for your computer
library(jagsUI)

# LOAD DATA #
harvest <- read.csv("AMWO harvest 2.csv",header=TRUE)
dim(harvest)
harvest[1,]

# read data into 51x2 matrices for jags (col 1 = Central, 2 = Eastern)
  # row 1 is 2013, row 51 is 1963 (reverse time)
HIP <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)
DSS <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)
stamps <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)

for (i in 1:51){
  HIP[i,2] <- harvest[i,13] 
  HIP[i,1] <- harvest[(i+51),13] 
  DSS[i,2] <- harvest[i,12] 
  DSS[i,1] <- harvest[(i+51),12] 
  stamps[i,2] <- harvest[i,22]
  stamps[i,1] <- harvest[(i+51),22]}

head(HIP)
tail(stamps)

# use mean and sd from first 10 years HIP (2004-2013) for year 1 priors
mean(HIP[1:10,1])
1/(sd(HIP[1:10,1]))^2 # precision
mean(HIP[1:10,2])
1/(sd(HIP[1:10,2]))^2

# use max duck stamps as lower bound on max hunters
max(stamps[,1],na.rm=TRUE)
max(stamps[,2],na.rm=TRUE)

########################################################
# Load WING and HARVEST data #
########################################################
harvest <- read.csv("AMWO harvest.csv",header=TRUE)
#NOTE: these data are missing 1997, 1998, 2014 and 2015!!

#bring in year
clean<-matrix(NA,nrow=102,ncol=1)
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

#this is the sum of all age-sex wing classes. 
clean[,9] <- rowSums(clean[,c(3:4,6:8)])    # NOTE: Adjust if don't want all unknown groups

colnames(clean)<-c("Year","Pop","AdF","AdM","AdU","JuvF","JuvM","JuvU","Total")

# Now need to make wings.age and wings.sex datasets
#clean[,10] <- rowSums(clean[,3:5])   # adult wings column; adjust categories included if needed
clean[,5] <- rowSums(clean[,6:8])   # juvenile wings column; adjust categories included if needed

colnames(clean)<-c("Year","Pop","AdF","AdM","JuvTot","JuvF","JuvM","JuvU","Total")

wings<-array(NA, dim=c(51,4,2))                # dims: t, ages (1=juv; 2=adult), p (1=eastern; 2=central)
wings[,1,1]<-clean$JuvTot[clean$Pop=="E"]      #juvenile eastern
wings[,2,1]<-clean$AdM[clean$Pop=="E"]         #adult male eastern
wings[,3,1]<-clean$AdF[clean$Pop=="E"]         #adult female eastern
wings[,4,1]<-clean$Total[clean$Pop=="E"]       # total wings eastern
wings[,1,2]<-clean$JuvTot[clean$Pop=="C"]      #juvenile central
wings[,2,2]<-clean$AdM[clean$Pop=="C"]         # adult male central
wings[,3,2]<-clean$AdF[clean$Pop=="C"]         # adult female central
wings[,4,2]<-clean$Total[clean$Pop=="C"]       # total wings central

# replace NAs in 1997 and 1998
wings[35:36,1,1] <- round(mean(wings[,1,1], na.rm=TRUE))
wings[35:36,2,1] <- round(mean(wings[,2,1], na.rm=TRUE))
wings[35:36,3,1] <- round(mean(wings[,3,1], na.rm=TRUE))
wings[35:36,4,1] <- round(mean(wings[,4,1], na.rm=TRUE))
wings[35:36,1,2] <- round(mean(wings[,1,2], na.rm=TRUE))
wings[35:36,2,2] <- round(mean(wings[,2,2], na.rm=TRUE))
wings[35:36,3,2] <- round(mean(wings[,3,2], na.rm=TRUE))
wings[35:36,4,2] <- round(mean(wings[,4,2], na.rm=TRUE))

#############################################################
# specify jags model to predict total annual harvest
#############################################################
sink("AMWO.harvest.wings.jags")
cat("
    model {
  # [,1] is central pop, [,2] is eastern population
  # prior for year 1 Harvest (reverse time, [1,]=2013)
    harvest[1] ~ dnorm(436095,1.822e-10) # mean, 1/sd2 from 2004-2013 data
    harvest[2] ~ dnorm(161459,9.551e-10)

  # round harvest estimates for first years
    Harvest[1,1] <- round(harvest[1])
    Harvest[1,2] <- round(harvest[2])

  # prior for max hunters in duck stamp years (up to 3-fold higher than max seen)
    max.hunters[1] ~ dunif(1000000,3000000)
    max.hunters[2] ~ dunif(500000,1500000)

  # priors for annual variation in true harvest (sigma.H) 
    # and estimated harvest (sd.HIP, sd.DSS)
    for (p in 1:2){
      sigma.H[p] ~ dunif(1,200000)
      tau.H[p] <- pow(sigma.H[p],-2)
      sd.HIP[p] ~ dunif(1,200000)
      prec.HIP[p] <- pow(sd.HIP[p],-2)
      sd.DSS[p] ~ dunif(1,200000)
      prec.DSS[p] <- pow(sd.DSS[p],-2)
    }

  # State process (estimate true Harvest in years 2-51)
    for (p in 1:2){
      for (i in 2:yrs){
      eps[i,p] ~ dnorm(0,tau.H[p])
      Harvest[i,p] <- trunc(Harvest[i-1,p] + eps[i,p]) # autoregressive model
      }
    }

# observation process, recoveries and harvest data
    for (p in 1:2){
      for (i in 1:15){
        HIP[i,p] ~ dnorm(Harvest[i,p],prec.HIP[p])} # end HIP years
      for (j in 13:50){
        frac[j,p] <- stamps[j,p]/max.hunters[p]
        DSS[j,p] ~ dnorm(frac[j,p]*Harvest[j,p],prec.DSS[p])} # end DSS yrs
    } # end p loop

    # Determining age-sex class proportions (pi's)

    # Summarize known wing samples by age and sex
    for (p in 1:2){

    pi[1,p] ~ dunif(0.1, 0.8)
    pi[2,p] ~ dunif(0.1, 0.8)
    pi[3,p] <- 1-(pi[1,p]+pi[2,p])

    for (t in 1:yrs){                                               # populations (1 eastern, 2 central)
    wings.class[t,,p] ~ dmulti(pi[,p], wings.total[t,p])          # pi is 3 proportions: 1=juvs, 2=males, 3=females    

    H[t,1,p] ~ dbin(pi[1,p], Harvest[rev[t],p])
    H[t,2,p] ~ dbin(pi[2,p], Harvest[rev[t],p])
    H[t,3,p] ~ dbin(pi[3,p], Harvest[rev[t],p])
    } #t
    } #p
  } # end jags model
    ",fill = TRUE)
sink()

# Bundle data  , add 1 to harvest, banding and recovery totals, per Seber's small-sample modification
wings.total <- wings[,4,]
wings.class <- wings[,1:3,]

#Reverse time vector
rev <- seq(51,1,-1)

bugs.data <- list(yrs=51, HIP=HIP, DSS=DSS, stamps=stamps, wings.total=wings.total, wings.class=wings.class,
                  rev=rev)

# Initial values (not all priors listed)

pi.inits <- matrix(nrow=3, ncol=2)
for (c in 1:3){
  for (p in 1:2){
    pi.inits[1,p] <- 0.46
    pi.inits[2,p] <- 0.21
    pi.inits[3,p] <- NA
  }}

inits <- function(){list(max.hunters=c(2200000,1000000), pi=pi.inits, harvest=runif(2,300000,500000))}  

# Parameters monitored
parameters <- c("Harvest","max.hunters","sd.DSS","sd.HIP","sigma.H","H","pi")

# MCMC settings
ni <- 2500000
nt <- 100        
nb <- 500000
nc <- 3

# call jags
AMWO.harvest.jags <- jagsUI(bugs.data, inits, parameters, "AMWO.harvest.wings.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)

print(AMWO.harvest.jags,digits=1)
