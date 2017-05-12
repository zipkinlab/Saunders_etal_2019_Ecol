#--------------------------------------------------------------------------------------------------------------
# Estimating annual harvest of AMWO during duck stamp years using annual sales as a covariate
#--------------------------------------------------------------------------------------------------------------

#rm(list=ls(all=TRUE)) # clear memory

# load R2WinBUGS, set working directory for data and output, and identify WinBUGS location for your computer
library(jagsUI)

# LOAD DATA #
harvest <- read.csv("AMWO harvest 2.csv",header=TRUE)
dim(harvest)
harvest[1,]

# read data into 51x2 matrices for jags (col 1 = Eastern, 2 = Central)
# row 1 is 2013, row 51 is 1963 (reverse time)
HIP <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)
DSS <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)
#stamps <- matrix(NA,nrow=51,ncol=2,byrow=TRUE)

for (i in 1:51){
  HIP[i,2] <- harvest[i,13] 
  HIP[i,1] <- harvest[(i+51),13] 
  DSS[i,2] <- harvest[i,12] 
  DSS[i,1] <- harvest[(i+51),12]}

#head(HIP)

##############################
# scale harvest estimates by dividing by 10000
##############################
HIP <- HIP/10000
DSS <- DSS/10000

# use mean and sd from first 10 years HIP (2004-2013) for year 1 priors. NOTE: Scaled now
mean(HIP[1:10,1])  #eastern
1/(sd(HIP[1:10,1]))^2 # precision
mean(HIP[1:10,2])  #central
1/(sd(HIP[1:10,2]))^2

# use max duck stamps as lower bound on max hunters
#max(stamps[,1],na.rm=TRUE)  #eastern
#max(stamps[,2],na.rm=TRUE)  #central

########################################################
# Load WING and HARVEST data #
########################################################
harvest <- read.csv("AMWO harvest.csv",header=TRUE)

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

wings<-array(NA, dim=c(51,4,2))                # dims: t, ages (1=juv; 2=adult male; 3=adult female; 4=total wings), p (1=eastern; 2=central)
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
    
    # [,1] is eastern pop, [,2] is central population
    # prior for year 1 Harvest (reverse time, [1,]=2013)
    harvest[1] ~ dnorm(16,9.551e-2) # mean, 1/sd2 from 2004-2013 data, eastern SCALED
    harvest[2] ~ dnorm(44,1.822e-2) # central SCALED
    
    # round harvest estimates for first years
    Harvest[1,1] <- round(harvest[1])
    Harvest[1,2] <- round(harvest[2])
    
    # priors for annual variation in true harvest (sigma.H) 
    # and estimated harvest (sd.HIP, sd.DSS)
    for (p in 1:2){
    sigma.H[p] ~ dunif(1,20)  # changed from 200000 for scaling
    tau.H[p] <- pow(sigma.H[p],-2)
    sd.HIP[p] ~ dunif(1,20)  # changed from 200000 for scaling
    prec.HIP[p] <- pow(sd.HIP[p],-2)
    sd.DSS[p] ~ dunif(1,20)  # changed from 200000 for scaling
    prec.DSS[p] <- pow(sd.DSS[p],-2)
    } #p
    
    # State process (estimate true Harvest in years 2-51)
    for (p in 1:2){
    for (i in 2:yrs){
    eps[i,p] ~ dnorm(0,tau.H[p])
    Harvest[i,p] <- trunc(Harvest[i-1,p] + eps[i,p]) # autoregressive model
    } #i
    } #p
    
    # observation process, recoveries and harvest data
    for (p in 1:2){
    for (i in 1:15){
    HIP[i,p] ~ dnorm(Harvest[i,p],prec.HIP[p])} #i  end HIP years
    for (j in 13:50){
    DSS[j,p] ~ dnorm(Harvest[j,p]/cor[p],prec.DSS[p])} #j end DSS yrs
    } # end p loop
    
    # Determining age-sex class proportions (pi's)
    for (p in 1:2){
    for (t in 1:yrs){
    pi[t,1,p] ~ dunif(0.1, 0.8)
    pi[t,2,p] ~ dunif(0.1, 0.8)
    pi[t,3,p] <- 1-(pi[t,1,p]+pi[t,2,p])

    # Summarize known wing samples by age and sex
    # populations (1 eastern, 2 central)
    wings.class[t,,p] ~ dmulti(pi[t,,p], wings.total[t,p])          # pi is 3 proportions: 1=juvs, 2=males, 3=females    
    H[t,1,p] <- pi[t,1,p]*Harvest[rev[t],p]
    H[t,2,p] <- pi[t,2,p]*Harvest[rev[t],p]
    H[t,3,p] <- pi[t,3,p]*Harvest[rev[t],p]
    } #t
    } #p
    } # end jags model
    ",fill = TRUE)
sink()

# Bundle data
wings.total <- wings[,4,]
wings.class <- wings[,1:3,]

#Reverse time vector
rev <- seq(51,1,-1)
cor <- c(3.1, 4.1)
bugs.data <- list(yrs=51, HIP=HIP, DSS=DSS, wings.total=wings.total, wings.class=wings.class,
                  rev=rev, cor=cor)

# Initial values (not all priors listed)

pi.inits <- array(dim=c(51,3,2))
for (c in 1:3){
  for (p in 1:2){
    for (t in 1:51){
    pi.inits[t,1,p] <- 0.46
    pi.inits[t,2,p] <- 0.21
    pi.inits[t,3,p] <- NA
  }}}

inits <- function(){list(pi=pi.inits, harvest=c(13,37))}   #harvest inits based on priors (means of 2004-2013), scaled from 130000 and 370000

# Parameters monitored
parameters <- c("Harvest","sd.DSS","sd.HIP","sigma.H","H","pi")  

# MCMC settings
ni <- 200000
nt <- 5          
nb <- 150000 
nc <- 3

# call jags
AMWO.harvest.wings.scaled.HPCC <- jagsUI(bugs.data, inits=inits, parameters, "AMWO.harvest.wings.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)
save(AMWO.harvest.wings.scaled.HPCC, file="AMWO_harvest_wings_scaled_model_HPCC.Rda")
#print(AMWO.harvest.wings.scaled.HPCC, digits=3)
#load(file="AMWO_harvest_wings_deterministic_model_HPCC.Rda")
#save(AMWO.harvest.wings.34cor, file="AMWO harvest wings deterministic practice.Rda")
