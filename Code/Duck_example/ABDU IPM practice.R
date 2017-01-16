# Integrated Population Model for American Black Ducks
# Todd Arnold, Univ. of Minnesota, arnol065@umn.edu (code sent 1/13/17)
rm(list=ls(all=TRUE)) # clear memory

# load R2WinBUGS, set working directory for data and output, and identify WinBUGS location for your computer
library(R2WinBUGS)  
library(R2jags)
library(jagsUI)
setwd("C:/Users/arnol065/Documents/BUGS/ABDU IPM") 
bugs.dir <- "C:/WinBUGS/WinBUGS14" 

# LOAD DATA #
harvest <- read.csv("ABDU harvest and wing data.csv",header=TRUE)
head(harvest)

# partition data sets by US and Canada (proportionate harvest differs by age and sex)
US.harvest <- matrix(0,nrow=20,ncol=13,byrow=TRUE)
Can.harvest <- matrix(0,nrow=20,ncol=13,byrow=TRUE)
US.wings <- matrix(0,nrow=20,ncol=8,byrow=TRUE)
Can.wings <- matrix(0,nrow=20,ncol=8,byrow=TRUE)
for (i in 1:20){
  for (k in 1:13){
    US.harvest[i,k] <- harvest[i,k]
    Can.harvest[i,k] <- harvest[i+20,k]
  }
  for (j in 1:6){
    US.wings[i,j] <- harvest[i,j+13]
    Can.wings[i,j] <- harvest[i+20,j+13]
  } #j
  US.wings[i,7] <- sum(US.wings[i,1:3]) # summarize known adults (including unknown sex)
  US.wings[i,8] <- sum(US.wings[i,4:6]) # summarize known juvs (including unknown sex)
  Can.wings[i,7] <- sum(Can.wings[i,1:3]) # summarize known adults (including unknown sex)
  Can.wings[i,8] <- sum(Can.wings[i,4:6]) # summarize known juvs (including unknown sex)
} #i

# banding data consists of 6 40x21 m-arrays, JuvF, AdF, UnkF (postseason bandings), JuvM, AdM, UnkM
# doubling of rows due to postseason (odd) and preseason (even) bandings
# no unk in preseason samples, so unk arrays include empty rows 
# m-arrays start in 1969 postseason (Jan-Mar) and move in 6 month increments to preseason 1988
# column 21 is banded birds never recovered (total releases is SUM ROW TOTAL) 

banding <- read.csv("ABDU shot only.csv",header=FALSE)
head(banding)

# create matrix of total bandings for Lincoln estimates
banded <- matrix(0,nrow=40,ncol=8,byrow=TRUE) # juv, ad, unk, sum, total bandings for each cohort in each banding period

# read data into individual m-arrays based on sex (F,M) and age at banding (juv,ad,unk)
marr.juvF.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.juvF.shot[i,j] <- banding[i,j]}
  banded[i,1]<-sum(marr.juvF.shot[i,])}

marr.adF.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.adF.shot[i,j] <- banding[(i+40),j]}
  banded[i,2]<-sum(marr.adF.shot[i,])}

marr.unkF.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.unkF.shot[i,j] <- banding[(i+80),j]}
  banded[i,3]<-sum(marr.unkF.shot[i,])
  banded[i,4]<-sum(banded[i,1:3])}  # total bandings

# merge adult and unknown for models with no age variation during summer S
marr.ad_unkF.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.ad_unkF.shot[i,j] <- marr.adF.shot[i,j]+marr.unkF.shot[i,j]}}

# regroup postseason bandings for Lincoln estimates
marr.combF.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE) 
for (i in 1:40){
  for (j in 1:21){
    marr.combF.shot[i,j] <- marr.juvF.shot[i,j]+marr.adF.shot[i,j]+marr.unkF.shot[i,j]
  }}

marr.juvM.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.juvM.shot[i,j] <- banding[(i+120),j]}
  banded[i,5]<-sum(marr.juvM.shot[i,])}

marr.adM.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.adM.shot[i,j] <- banding[(i+160),j]}
  banded[i,6]<-sum(marr.adM.shot[i,])}

marr.unkM.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.unkM.shot[i,j] <- banding[(i+200),j]}
  banded[i,7]<-sum(marr.unkM.shot[i,])
  banded[i,8]<-sum(banded[i,5:7])}

marr.ad_unkM.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE)
for (i in 1:40){
  for (j in 1:21){
    marr.ad_unkM.shot[i,j] <- marr.adM.shot[i,j]+marr.unkM.shot[i,j]}}

marr.combM.shot <- matrix(0,nrow=40,ncol=21,byrow=TRUE) # regroup postseason bandings for Lincoln estimates
for (i in 1:40){
  for (j in 1:21){
    marr.combM.shot[i,j] <- marr.juvM.shot[i,j]+marr.adM.shot[i,j]+marr.unkM.shot[i,j]
  }}

# inspect banding totals
tail(banded)

#----------------------------------------------------------------------------------------
# Part 1. 
# Generate Lincoln estimates for each age sex class to use as population estimates (y[t])
#----------------------------------------------------------------------------------------

# create matrix to store Lincoln estimates and SD as output for part 2
Lincoln.out <- matrix(NA, nrow = 20, ncol = 18) 

# specify bugs model
sink("Lincoln.bug")
cat("
    model {
    
    # Priors and constraints for true proportions of population in each age-sex class 
    # given constraint of summing to 1, these 3 proportions fully describe the data each year
    pi.US.juv.mu ~ dunif(-1, 1)       # prior for mean proportion of juvs in US Harvest
    pi.US.juvF.mu ~ dunif(-1, 1)      # mean proportion of Females in juvenile US Harvest
    pi.US.adF.mu ~ dunif(-1, 1)       # mean proportion of Females in the adult US harvest
    pi.Can.juv.mu ~ dunif(-1, 3)       
    pi.Can.juvF.mu ~ dunif(-1, 1)     
    pi.Can.adF.mu ~ dunif(-1, 1)      
    
    # Priors and constraints, for logit link on f
    f.pre.juvF.mu ~ dunif(-6, 0)       # priors for mean Brownie recovery (bounded 0, 0.5)
    f.pre.adF.mu ~ dunif(-6, 0)        
    f.post.combF.mu ~ dunif(-6, 0)     # harvest is mixed, so pooled recov rate  
    f.pre.juvM.mu ~ dunif(-6, 0)       
    f.pre.adM.mu ~ dunif(-6, 0)        
    f.post.combM.mu ~ dunif(-6, 0) 
    
    report <- 0.506                 # reporting rate (Conroy & Blandin, fixed, can estimate in joint analysis)
    harv.adj <- 0.72                 # harvest adjustment (Padding & Royle, adjust for positive bias in US harvest)
    
    pi.US.juv.sd ~ dunif(0, 1)       # prior for annual SD
    pi.US.juvF.sd ~ dunif(0, 1)      
    pi.US.adF.sd ~ dunif(0, 1)       
    pi.Can.juv.sd ~ dunif(0, 1)       
    pi.Can.juvF.sd ~ dunif(0, 1)      
    pi.Can.adF.sd ~ dunif(0, 1)       
    
    # priors for annual SD of annual recovery rate (logit scale)
    f.pre.juvF.sd ~ dunif(0,2)       
    f.pre.adF.sd ~ dunif(0,2)        
    f.post.combF.sd ~ dunif(0,2)       
    f.pre.juvM.sd ~ dunif(0,2)       
    f.pre.adM.sd ~ dunif(0,2)        
    f.post.combM.sd ~ dunif(0,2)       
    
    pi.US.juv.tau <- pow(pi.US.juv.sd,-2)       # convert to tau
    pi.US.juvF.tau <- pow(pi.US.juvF.sd,-2)      
    pi.US.adF.tau <- pow(pi.US.adF.sd,-2)       
    pi.Can.juv.tau <- pow(pi.Can.juv.sd,-2)       
    pi.Can.juvF.tau <- pow(pi.Can.juvF.sd,-2)      
    pi.Can.adF.tau <- pow(pi.Can.adF.sd,-2)       
    # express variances as precision (1/SD^2)
    f.pre.juvF.tau <- pow(f.pre.juvF.sd,-2)
    f.pre.adF.tau  <- pow(f.pre.adF.sd,-2)
    f.pre.juvM.tau  <- pow(f.pre.juvM.sd,-2)
    f.pre.adM.tau  <- pow(f.pre.adM.sd,-2)
    f.post.combF.tau  <- pow(f.post.combF.sd,-2)
    f.post.combM.tau  <- pow(f.post.combM.sd,-2)
    
    # Generate annual parameter estimates 
    for (y in 1:yrs){
    logit.pi.US.juv[y] ~ dnorm(pi.US.juv.mu,pi.US.juv.tau)
    logit.pi.US.juvF[y] ~ dnorm(pi.US.juvF.mu,pi.US.juvF.tau)
    logit.pi.US.adF[y] ~ dnorm(pi.US.adF.mu,pi.US.adF.tau)
    logit.pi.Can.juv[y] ~ dnorm(pi.Can.juv.mu,pi.Can.juv.tau)
    logit.pi.Can.juvF[y] ~ dnorm(pi.Can.juvF.mu,pi.Can.juvF.tau)
    logit.pi.Can.adF[y] ~ dnorm(pi.Can.adF.mu,pi.Can.adF.tau)
    logit.f.pre.juvF[y] ~ dnorm(f.pre.juvF.mu,f.pre.juvF.tau)
    logit.f.pre.adF[y] ~ dnorm(f.pre.adF.mu,f.pre.adF.tau)
    logit.f.pre.juvM[y] ~ dnorm(f.pre.juvM.mu,f.pre.juvM.tau)
    logit.f.pre.adM[y] ~ dnorm(f.pre.adM.mu,f.pre.adM.tau)
    logit.f.post.combF[y] ~ dnorm(f.post.combF.mu,f.post.combF.tau)
    logit.f.post.combM[y] ~ dnorm(f.post.combM.mu,f.post.combM.tau)
    logit(pi.US.juv[y]) <- logit.pi.US.juv[y]
    logit(pi.US.juvF[y]) <- logit.pi.US.juvF[y]
    logit(pi.US.adF[y]) <- logit.pi.US.adF[y]
    logit(pi.Can.juv[y]) <- logit.pi.Can.juv[y]
    logit(pi.Can.juvF[y]) <- logit.pi.Can.juvF[y]
    logit(pi.Can.adF[y]) <- logit.pi.Can.adF[y]
    logit(f.pre.juvF[y]) <- logit.f.pre.juvF[y]
    logit(f.pre.adF[y]) <- logit.f.pre.adF[y]
    logit(f.pre.juvM[y]) <- logit.f.pre.juvM[y]
    logit(f.pre.adM[y]) <- logit.f.pre.adM[y]
    logit(f.post.combF[y]) <- logit.f.post.combF[y]
    logit(f.post.combM[y]) <- logit.f.post.combM[y]
    }
    
    # Summarize known wing samples
    for (y in 1:yrs){
    US.wings[y,8] ~ dbin(pi.US.juv[y],US.known.age[y])
    US.wings[y,4] ~ dbin(pi.US.juvF[y],US.sexed.juv[y])
    US.wings[y,1] ~ dbin(pi.US.adF[y],US.sexed.ad[y])
    Can.wings[y,8] ~ dbin(pi.Can.juv[y],Can.known.age[y])
    Can.wings[y,4] ~ dbin(pi.Can.juvF[y],Can.sexed.juv[y])
    Can.wings[y,1] ~ dbin(pi.Can.adF[y],Can.sexed.ad[y])
    }
    
    # Generate component estimates of harvest
    for (y in 1:yrs){
    Harv.JuvF.Can[y] <- pi.Can.juv[y]*pi.Can.juvF[y]*Can.Harvest[y,12]
    Harv.JuvM.Can[y] <- pi.Can.juv[y]*(1-pi.Can.juvF[y])*Can.Harvest[y,12]
    Harv.AdF.Can[y] <- (1-pi.Can.juv[y])*pi.Can.adF[y]*Can.Harvest[y,12]
    Harv.AdM.Can[y] <- (1-pi.Can.juv[y])*(1-pi.Can.adF[y])*Can.Harvest[y,12]
    Harv.JuvF.US[y] <- pi.US.juv[y]*pi.US.juvF[y]*US.Harvest[y,12]*harv.adj
    Harv.JuvM.US[y] <- pi.US.juv[y]*(1-pi.US.juvF[y])*US.Harvest[y,12]*harv.adj
    Harv.AdF.US[y] <- (1-pi.US.juv[y])*pi.US.adF[y]*US.Harvest[y,12]*harv.adj
    Harv.AdM.US[y] <- (1-pi.US.juv[y])*(1-pi.US.adF[y])*US.Harvest[y,12]*harv.adj
    Harv.JuvF.Total[y] <- Harv.JuvF.Can[y]+Harv.JuvF.US[y]
    Harv.AdF.Total[y] <- Harv.AdF.Can[y]+Harv.AdF.US[y]
    Harv.JuvM.Total[y] <- Harv.JuvM.Can[y]+Harv.JuvM.US[y]
    Harv.AdM.Total[y] <- Harv.AdM.Can[y]+Harv.AdM.US[y]
    }
    
    # Generate direct recovery rates
    for (y in 1:yrs){
    marr.juvF.shot[y*2,y] ~ dbin(f.pre.juvF[y],banded[y*2,1])
    marr.adF.shot[y*2,y] ~ dbin(f.pre.adF[y],banded[y*2,2])
    marr.combF.shot[y*2-1,y] ~ dbin(f.post.combF[y],banded[y*2-1,4])
    marr.juvM.shot[y*2,y] ~ dbin(f.pre.juvM[y],banded[y*2,5])
    marr.adM.shot[y*2,y] ~ dbin(f.pre.adM[y],banded[y*2,6])
    marr.combM.shot[y*2-1,y] ~ dbin(f.post.combM[y],banded[y*2-1,8])
    }
    
    # Generate Lincoln estimates and population structure as derived parameters 
    for (y in 1:yrs){
    N.juvF.fall[y] <- (report*Harv.JuvF.Total[y])/f.pre.juvF[y]
    N.adF.fall[y] <- (report*Harv.AdF.Total[y])/f.pre.adF[y]
    N.combF.spring[y] <- (report*Harv.AdF.Total[y])/f.post.combF[y]
    N.juvM.fall[y] <- (report*Harv.JuvM.Total[y])/f.pre.juvM[y]
    N.adM.fall[y] <- (report*Harv.AdM.Total[y])/f.pre.adM[y]
    N.combM.spring[y] <- (report*Harv.AdM.Total[y])/f.post.combM[y]
    N.total.fall[y] <- N.juvF.fall[y] + N.juvM.fall[y] + N.adF.fall[y] + N.adM.fall[y]
    N.total.spring[y] <- N.combF.spring[y] + N.combM.spring[y]
    fecundity[y] <- 0.5*(N.juvF.fall[y]+N.juvM.fall[y])/N.adF.fall[y]
    sex.ratio.juv[y] <- N.juvM.fall[y]/(N.juvF.fall[y]+N.juvM.fall[y])
    sex.ratio.ad.fall[y] <- N.adM.fall[y]/(N.adF.fall[y]+N.adM.fall[y])
    sex.ratio.comb.spring[y] <- N.combM.spring[y]/(N.combF.spring[y]+N.combM.spring[y])
    
    }
    
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data  # marr.unkF.shot=marr.unkF.shot, marr.unkM.shot=marr.unkM.shot, 
bugs.data <- list(US.Harvest=US.harvest, Can.Harvest=Can.harvest, US.wings=US.wings, Can.wings=Can.wings, banded=banded, 
                  marr.juvF.shot=marr.juvF.shot, marr.adF.shot=marr.adF.shot, marr.combF.shot=marr.combF.shot, 
                  marr.juvM.shot=marr.juvM.shot, marr.adM.shot=marr.adM.shot, marr.combM.shot=marr.combM.shot, yrs = 20,
                  US.known.age = US.wings[,7]+US.wings[,8], US.sexed.juv = US.wings[,4]+US.wings[,5],
                  US.sexed.ad = US.wings[,1]+US.wings[,2], Can.known.age = Can.wings[,7]+Can.wings[,8],
                  Can.sexed.juv = Can.wings[,4]+Can.wings[,5], Can.sexed.ad = Can.wings[,1]+Can.wings[,2])


# Initial values (not all priors listed)
inits <- function(){list(pi.US.juv.mu = runif(1, -1, 1), pi.US.juvF.mu = runif(1, -1, 1), pi.US.adF.mu = runif(1, -1, 1),
                         pi.US.juv.sd = runif(1, 0, 1), pi.US.juvF.sd = runif(1, 0, 1), pi.US.adF.sd = runif(1,0,1), 
                         pi.Can.juv.mu = runif(1, -1, 2), pi.Can.juvF.mu = runif(1, -1, 1), pi.Can.adF.mu = runif(1, -1, 1),
                         pi.Can.juv.sd = runif(1, 0, 1), pi.Can.juvF.sd = runif(1, 0, 1), pi.Can.adF.sd = runif(1,0,1))}  

# Parameters monitored
parameters <- c("Harv.JuvF.Total","Harv.JuvM.Total","Harv.AdF.Total","Harv.AdM.Total",
                "N.juvF.fall","N.juvM.fall","N.adF.fall","N.adM.fall","N.total.fall","N.combF.spring","N.combM.spring","N.total.spring",
                "fecundity","sex.ratio.juv","sex.ratio.ad.fall","sex.ratio.comb.spring")

# MCMC settings (90 sec, converges rapidly, great mixing)
ni <- 25000
nt <- 6        #10,000 posterior samples
nb <- 5000
nc <- 3

# Call WinBUGS from R (82 sec, this part is very fast, nice mixing, converges well)
ABDU.Lincoln.bug <- bugs(bugs.data, inits, parameters, "Lincoln.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
print(ABDU.Lincoln.bug, digits = 4)

ABDU.Lincoln.jags <- jagsUI(bugs.data, inits, parameters, "Lincoln.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(ABDU.Lincoln.jags, digits = 3)

names(ABDU.Lincoln.bug$sims.list)
hist(ABDU.Lincoln.bug$sims.list$N.juvF.fall[,1])
hist(ABDU.Lincoln.bug$sims.list$sex.ratio.ad.fall[,19])
hist(ABDU.Lincoln.bug$sims.list$fecundity[,19])


# save Lincoln output to data file (might be easier if I truncated N and coverted SD to tau here)
for (y in 1:20){
  Lincoln.out[y,1]=trunc(ABDU.Lincoln.bug$mean$N.juvF.fall[y])
  Lincoln.out[y,2]=trunc(ABDU.Lincoln.bug$sd$N.juvF.fall[y])
  Lincoln.out[y,3]=trunc(ABDU.Lincoln.bug$mean$N.juvM.fall[y])
  Lincoln.out[y,4]=trunc(ABDU.Lincoln.bug$sd$N.juvM.fall[y])
  Lincoln.out[y,5]=trunc(ABDU.Lincoln.bug$mean$N.adF.fall[y])
  Lincoln.out[y,6]=trunc(ABDU.Lincoln.bug$sd$N.adF.fall[y])
  Lincoln.out[y,7]=trunc(ABDU.Lincoln.bug$mean$N.adM.fall[y])
  Lincoln.out[y,8]=trunc(ABDU.Lincoln.bug$sd$N.adM.fall[y])
  Lincoln.out[y,9]=trunc(ABDU.Lincoln.bug$mean$N.combF.spring[y])
  Lincoln.out[y,10]=trunc(ABDU.Lincoln.bug$sd$N.combF.spring[y])
  Lincoln.out[y,11]=trunc(ABDU.Lincoln.bug$mean$N.combM.spring[y])
  Lincoln.out[y,12]=trunc(ABDU.Lincoln.bug$sd$N.combM.spring[y])
  Lincoln.out[y,13]=ABDU.Lincoln.bug$mean$fecundity[y]
  Lincoln.out[y,14]=ABDU.Lincoln.bug$sd$fecundity[y]
  Lincoln.out[y,15]=trunc(ABDU.Lincoln.bug$mean$N.adF.fall[y]+ABDU.Lincoln.bug$mean$N.adF.fall[y])
  Lincoln.out[y,16]=sqrt((ABDU.Lincoln.bug$sd$N.adF.fall[y])^2+(ABDU.Lincoln.bug$sd$N.adF.fall[y])^2)
  Lincoln.out[y,17]=trunc(ABDU.Lincoln.bug$mean$N.adM.fall[y]+ABDU.Lincoln.bug$mean$N.adM.fall[y])
  Lincoln.out[y,18]=sqrt((ABDU.Lincoln.bug$sd$N.adM.fall[y])^2+(ABDU.Lincoln.bug$sd$N.adM.fall[y])^2)
}

head(Lincoln.out)
write.csv(Lincoln.out,file="Lincoln estimates.csv")

# summarize initial year log values for use as priors
log(Lincoln.out[1,9])  # prior for adult F spring pop
log(Lincoln.out[1,9])-log(Lincoln.out[1,9]-Lincoln.out[1,10]) # lower SD on log scale (lower is larger)
1/(log(Lincoln.out[1,9])-log(Lincoln.out[1,9]-Lincoln.out[1,10]))^2 # precision
log(Lincoln.out[1,11]) # prior for adult M spring pop
log(Lincoln.out[1,11])-log(Lincoln.out[1,11]-Lincoln.out[1,12]) #lower SD on log scale (lower is larger)
1/(log(Lincoln.out[1,11])-log(Lincoln.out[1,11]-Lincoln.out[1,12]))^2 # precision

#----------------------------------------------------------------------------------------
# Part 1.b 
# Generate fecundity estimates from age-ratios at capture
# Note that these are biased unrealistically high, don't use them as estimates of F
#----------------------------------------------------------------------------------------
live.recaps <- read.csv("live.recaptures.csv",header=TRUE)
head(live.recaps)

# Specify model in BUGS language
sink("age_ratio.bug")
cat("
    model {
    
    # Priors and constraints
    p.jf ~ dunif(0, 1)         # Vague uniform prior for juv. first year recapture rate
    p.af ~ dunif(0, 1)         
    p.jm ~ dunif(0, 1)         
    p.am ~ dunif(0, 1)         
    
    # summarize total bandings
    for (y in 1:yrs){
    jf.banded[y] <- banded[y*2,1]
    af.banded[y] <- banded[y*2,2]
    jm.banded[y] <- banded[y*2,5]
    am.banded[y] <- banded[y*2,6]
    }
    
    #N.banded.jf <- sum(jf.banded[])
    #N.banded.af <- sum(af.banded[])
    #N.banded.jm <- sum(jm.banded[])
    #N.banded.am <- sum(am.banded[])
    
    # vulnerability adjustments, single value for all years
    recap.jf ~ dbin(p.jf,N.banded.jf)
    recap.af ~ dbin(p.af,N.banded.af)
    recap.jm ~ dbin(p.jm,N.banded.jm)
    recap.am ~ dbin(p.am,N.banded.am)
    vuln.jf <- p.jf/p.af
    vuln.jm <- p.jm/p.af
    
    # Annual fecundity estimates
    for (y in 1:yrs){
    F.jf[y] <- (jf.banded[y]/af.banded[y])/(p.jf/p.af)
    F.jm[y] <- (jm.banded[y]/af.banded[y])/(p.jm/p.af)
    } #y
    }
    ",fill = TRUE)
sink()

# Bundle data  # marr.unkF.shot=marr.unkF.shot, marr.unkM.shot=marr.unkM.shot, 
bugs.data <- list(banded=banded, recap.jf=sum(live.recaps[,2]), recap.af=sum(live.recaps[,3]), 
                  recap.jm=sum(live.recaps[,4]), recap.am=sum(live.recaps[,5]), yrs = 20,
                  N.banded.jf = sum(banded[,1]),
                  N.banded.af = sum(banded[,2]),
                  N.banded.jm = sum(banded[,5]),
                  N.banded.am = sum(banded[,6]))

# Initial values 
inits <- function(){list(p.jf = runif(1, 0, 1), p.af = runif(1, 0, 1), 
                         p.jm = runif(1, 0, 1), p.am = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("p.jf","p.af","p.jm","p.am","vuln.jf","vuln.jm","F.jf","F.jm")

# MCMC settings (90 sec, converges rapidly, good mixing)
ni <- 25000
nt <- 1
nb <- 5000
nc <- 3

# Call WinBUGS from R (82 sec, this part is very fast, nice mixing, converges well, but jags version doesn't work)
age_ratio.bug <- bugs(bugs.data, inits, parameters, "age_ratio.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
print(age_ratio.bug, digits = 4)

age_ratio.jags <- jagsUI(bugs.data, inits, parameters, "age_ratio.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(age_ratio.jags, digits = 4)


age_ratio.out <- matrix(NA, nrow = 20, ncol = 2) 
for (y in 1:20){
  age.ratio.out[y,1]=age_ratio.bug$mean$F.jf[y]
  age.ratio.out[y,2]=ABDU.Lincoln.bug$sd$F.jf[y]
  age.ratio.out[y,1]=age_ratio.bug$mean$F.jm[y]
  age.ratio.out[y,2]=ABDU.Lincoln.bug$sd$F.jm[y]
}  

#################################################
# Part 2: simple state-space IPM on Lincoln data
# log-normal version
################################################

sink("ABDU.log_SS.bug")  
cat("
    model {
    
    N.adF.spring[1] ~ dnorm(13.77,167) # observed initial Lincoln estimates and precision (log scale)
    N.adM.spring[1] ~ dnorm(13.85,286)
    
    m.adF.sum.mu ~ dunif(-1,0)    # instantaneous mortality rate, adF sum to fall 
    m.adM.sum.mu ~ dunif(-1,0)
    m.adF.win.mu ~ dunif(-1,0)
    m.adM.win.mu ~ dunif(-1,0)
    f.adF_juvF.mu ~ dunif(-2,2)      # fecundity rate, AdF to JuvF
    f.adF_juvM.mu ~ dunif(-2,2)      # fecundity rate, AdF to JuvM
    
    m.adF.sum.sd ~ dunif(0,0.5)     # sd of vital rates on log scale 
    m.adM.sum.sd ~ dunif(0,0.5)
    m.adF.win.sd ~ dunif(0,0.5)
    m.adM.win.sd ~ dunif(0,0.5)
    f.adF_juvF.sd ~ dunif(0,0.5)    
    f.adF_juvM.sd ~ dunif(0,0.5)    
    
    m.adF.sum.tau <- pow(m.adF.sum.sd,-2) 
    m.adM.sum.tau <- pow(m.adM.sum.sd,-2)
    m.adF.win.tau <- pow(m.adF.win.sd,-2)
    m.adM.win.tau <- pow(m.adM.win.sd,-2)
    f.adF_juvF.tau <- pow(f.adF_juvF.sd,-2) 
    f.adF_juvM.tau <- pow(f.adF_juvM.sd,-2)
    
    for (y in 1:yrs){
    m.adF.sum[y] ~ dnorm(m.adF.sum.mu,m.adF.sum.tau)
    m.adM.sum[y] ~ dnorm(m.adM.sum.mu,m.adM.sum.tau)
    m.adF.win[y] ~ dnorm(m.adF.win.mu,m.adF.win.tau) 
    m.adM.win[y] ~ dnorm(m.adM.win.mu,m.adM.win.tau) 
    f.adF_juvF[y] ~ dnorm(f.adF_juvF.mu,f.adF_juvF.tau) 
    f.adF_juvM[y] ~ dnorm(f.adF_juvM.mu,f.adF_juvM.tau) 
    }
    
    # Likelihood
    # State process (estimate fall parameters only in year 1)
    N.adF.fall[1] <- N.adF.spring[1] + m.adF.sum[1]
    N.adM.fall[1] <- N.adM.spring[1] + m.adM.sum[1]
    N.juvF[1] <- N.adF.fall[1] + f.adF_juvF[1]
    N.juvM[1] <- N.adF.fall[1] + f.adF_juvM[1]
    
    for (t in 2:yrs){
    N.adF.spring[t] <- log(exp(N.adF.fall[t-1]) + exp(N.juvF[t-1])) + m.adF.win[t-1]
    N.adM.spring[t] <- log(exp(N.adM.fall[t-1]) + exp(N.juvM[t-1])) + m.adM.win[t-1]
    N.adF.fall[t] <- N.adF.spring[t] + m.adF.sum[t]
    N.adM.fall[t] <- N.adM.spring[t] + m.adM.sum[t]
    N.juvF[t] <- N.adF.fall[t] + f.adF_juvF[t]
    N.juvM[t] <- N.adF.fall[t] + f.adF_juvM[t]
    }
    
    # Observation process
    for (t in 1:yrs) {
    juvF.obs[t] ~ dnorm(N.juvF[t], juvF.tau[t]) 
    juvM.obs[t] ~ dnorm(N.juvM[t], juvM.tau[t]) 
    adF.fall.obs[t] ~ dnorm(N.adF.fall[t], adF.fall.tau[t]) 
    adM.fall.obs[t] ~ dnorm(N.adM.fall[t], adM.fall.tau[t]) 
    adF.spring.obs[t] ~ dnorm(N.adF.spring[t], adF.spring.tau[t]) 
    adM.spring.obs[t] ~ dnorm(N.adM.spring[t], adM.spring.tau[t]) 
    }
    
    # back transform derived parameters
    for (t in 1:yrs) {
    juvF.pred[t] <- exp(N.juvF[t]) 
    juvM.pred[t] <- exp(N.juvM[t]) 
    adF.fall.pred[t] <- exp(N.adF.fall[t]) 
    adM.fall.pred[t] <- exp(N.adM.fall[t]) 
    adF.spring.pred[t] <- exp(N.adF.spring[t]) 
    adM.spring.pred[t] <- exp(N.adM.spring[t]) 
    total.fall.pred[t] <- juvF.pred[t]+juvM.pred[t]+adF.fall.pred[t]+adM.fall.pred[t]
    total.spring.pred[t] <- adF.spring.pred[t]+adM.spring.pred[t]
    }
    
    for (t in 1:(yrs-1)) {
    lambda.spring[t] <- total.spring.pred[t+1]/total.spring.pred[t]
    lambda.fall[t] <- total.fall.pred[t+1]/total.fall.pred[t]
    }
    
    }    
    ",fill=TRUE)
sink()

# Bundle data 
bugs.data <- list(juvF.obs=log(Lincoln.out[,1]), juvM.obs=log(Lincoln.out[,3]), adF.fall.obs=log(Lincoln.out[,5]),adM.fall.obs=log(Lincoln.out[,7]),
                  adF.spring.obs=log(Lincoln.out[,9]),adM.spring.obs=log(Lincoln.out[,11]),
                  juvF.tau=1/(log(Lincoln.out[,1])-log(Lincoln.out[,1]-Lincoln.out[,2]))^2,
                  juvM.tau=1/(log(Lincoln.out[,3])-log(Lincoln.out[,3]-Lincoln.out[,4]))^2,
                  adF.fall.tau=1/(log(Lincoln.out[,5])-log(Lincoln.out[,5]-Lincoln.out[,6]))^2,
                  adM.fall.tau=1/(log(Lincoln.out[,7])-log(Lincoln.out[,7]-Lincoln.out[,8]))^2,
                  adF.spring.tau=1/(log(Lincoln.out[,9])-log(Lincoln.out[,9]-Lincoln.out[,10]))^2,
                  adM.spring.tau=1/(log(Lincoln.out[,11])-log(Lincoln.out[,11]-Lincoln.out[,12]))^2,yrs=20)

# Initial values  
inits <- function() {list(m.adF.sum.mu = runif(1,-1,0), m.adM.sum.mu = runif(1,-1,0), 
                          m.adF.win.mu = runif(1,-1,0), m.adM.win.mu = runif(1,-1,0),
                          f.adF_juvF.mu = runif(1,-0.5,1), f.adF_juvM.mu = runif(1,-0.5,1),  
                          m.adF.sum.sd = runif(1,0.25,0.5), m.adM.sum.sd = runif(1,0.25,0.5), 
                          m.adF.win.sd = runif(1,0.25,0.5), m.adM.win.sd = runif(1,0.25,0.5),
                          f.adF_juvF.sd = runif(1,0.25,0.5), f.adF_juvM.sd = runif(1,0.25,0.5),   
                          N.adF.spring=c(rnorm(1,13.77,0.077),rep(NA,19)),
                          N.adM.spring=c(rnorm(1,13.85,0.059),rep(NA,19)))}

# Parameters monitored   
parameters <- c("m.adF.sum.mu", "m.adF.win.mu","m.adM.sum.mu", "m.adM.win.mu","f.adF_juvF.mu","f.adF_juvM.mu", 
                "m.adF.sum.sd", "m.adF.win.sd","m.adM.sum.sd", "m.adM.win.sd","f.adF_juvF.sd","f.adF_juvM.sd", 
                "m.adF.sum","m.adF.win","m.adM.sum","m.adM.win","f.adF_juvF","f.adF_juvM",
                "juvF.pred","juvM.pred","adM.fall.pred","adM.spring.pred", "adF.fall.pred","adF.spring.pred",
                "total.fall.pred","total.spring.pred","lambda.spring","lambda.fall")

# MCMC settings (150000 w 50000 burn gives good convergence on all parms R-hat <= 0.01)
# token ni to keep run time < 1 min
ni <- 1200
nt <- 1
nb <- 200
nc <- 3

# Call WinBUGS from R (bugs version hangs and never initiates)
#ABDU.log_SS.out <- bugs(bugs.data, inits, parameters, "ABDU_SS.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,debug=TRUE,bugs.directory=bugs.dir,working.directory=getwd())
#print(ABDU.log_SS.out,digits=2)

# jags version works well
ABDU.log_SS.jags <- jagsUI(bugs.data, inits, parameters, "ABDU.log_SS.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
print(ABDU.log_SS.jags,digits=3)

##########################
# Part 2b: simple State Space IPM Model on Lincoln data
# Poisson version
##########################
print(Lincoln.out[1,9]/1000) # prior for adF spring
print(Lincoln.out[1,11]/1000) # prior for adM spring

sink("ABDU.Pois_SS.bug")  
cat("
    model {
    
    N.adF.spring[1] ~ dpois(958.3)    # observed initial Lincoln estimates/1000
    N.adM.spring[1] ~ dpois(1033.9)
    
    S.adF.sum.mu ~ dunif(0,4)    # priors for seasonal survival rates (logit scale)
    S.adM.sum.mu ~ dunif(0,4)
    S.adF.win.mu ~ dunif(0,4)
    S.adM.win.mu ~ dunif(0,4)
    fecun.juvF.mu ~ dunif(1,2)      # priors for fecundity rates (sex spec because more F than M)
    fecun.juvM.mu ~ dunif(1,2)      
    
    S.adF.sum.sd ~ dunif(0,2)     # sd of vital rates
    S.adM.sum.sd ~ dunif(0,2)
    S.adF.win.sd ~ dunif(0,2)
    S.adM.win.sd ~ dunif(0,2)
    fecun.juvF.sd ~ dunif(0,1)    
    fecun.juvM.sd ~ dunif(0,1)    
    
    S.adF.sum.tau <- pow(S.adF.sum.sd,-2) 
    S.adM.sum.tau <- pow(S.adM.sum.sd,-2)
    S.adF.win.tau <- pow(S.adF.win.sd,-2)
    S.adM.win.tau <- pow(S.adM.win.sd,-2)
    fecun.juvF.tau <- pow(fecun.juvF.sd,-2) 
    fecun.juvM.tau <- pow(fecun.juvM.sd,-2)
    
    for (y in 1:yrs){
    logit.S.adF.sum[y] ~ dnorm(S.adF.sum.mu,S.adF.sum.tau)
    logit.S.adM.sum[y] ~ dnorm(S.adM.sum.mu,S.adM.sum.tau)
    logit.S.adF.win[y] ~ dnorm(S.adF.win.mu,S.adF.win.tau) 
    logit.S.adM.win[y] ~ dnorm(S.adM.win.mu,S.adM.win.tau) 
    logit(S.adF.sum[y]) <- logit.S.adF.sum[y] 
    logit(S.adM.sum[y]) <- logit.S.adM.sum[y] 
    logit(S.adF.win[y]) <- logit.S.adF.win[y] 
    logit(S.adM.win[y]) <- logit.S.adM.win[y] 
    Fecundity.juvF[y] ~ dnorm(fecun.juvF.mu,fecun.juvF.tau) 
    fecun.juvF[y] <-max(0,Fecundity.juvF[y])                # avoiding I(0,) for jags
    Fecundity.juvM[y] ~ dnorm(fecun.juvM.mu,fecun.juvM.tau) 
    fecun.juvM[y] <-max(0,Fecundity.juvM[y])
    }
    
    # Likelihood
    # State process (estimate fall parameters only in year 1)
    N.adF.fall[1] <- N.adF.spring[1] * S.adF.sum[1]
    N.adM.fall[1] <- N.adM.spring[1] * S.adM.sum[1]
    N.juvF[1] <- N.adF.fall[1] * fecun.juvF[1]
    N.juvM[1] <- N.adF.fall[1] * fecun.juvM[1]
    
    for (t in 2:yrs){
    N.adF.spring[t] <- (N.adF.fall[t-1] + N.juvF[t-1]) * S.adF.win[t-1]
    N.adM.spring[t] <- (N.adM.fall[t-1] + N.juvM[t-1]) * S.adM.win[t-1]
    N.adF.fall[t] <- N.adF.spring[t] * S.adF.sum[t]
    N.adM.fall[t] <- N.adM.spring[t] * S.adM.sum[t]
    N.juvF[t] <- N.adF.fall[t] * fecun.juvF[t]
    N.juvM[t] <- N.adF.fall[t] * fecun.juvM[t]
    lambda.spring[t-1] <- N.adF.spring[t]/N.adF.spring[t-1]
    lambda.fall[t-1] <- N.adF.fall[t]/N.adF.fall[t-1]
    }
    
    # Observation process
    for (t in 1:yrs) {
    juvF.obs[t] ~ dpois(N.juvF[t]) 
    juvM.obs[t] ~ dpois(N.juvM[t]) 
    adF.fall.obs[t] ~ dpois(N.adF.fall[t]) 
    adM.fall.obs[t] ~ dpois(N.adM.fall[t]) 
    adF.spring.obs[t] ~ dpois(N.adF.spring[t]) 
    adM.spring.obs[t] ~ dpois(N.adM.spring[t]) 
    }
    }    
    ",fill=TRUE)
sink()

# Bundle data 
bugs.data <- list(juvF.obs=trunc(Lincoln.out[,1]/1000), juvM.obs=trunc(Lincoln.out[,3]/1000), adF.fall.obs=trunc(Lincoln.out[,5]/1000),
                  adM.fall.obs=trunc(Lincoln.out[,7]/1000), adF.spring.obs=trunc(Lincoln.out[,9]/1000),adM.spring.obs=trunc(Lincoln.out[,11]/1000), yrs=20)

# Initial values  
# 
inits <- function() {list(S.adF.sum.mu = runif(1,0,4), S.adM.sum.mu = runif(1,0,4), 
                          S.adF.win.mu = runif(1,0,4), S.adM.win.mu = runif(1,0,4),
                          fecun.juvF.mu = runif(1,1.4,1.6), fecun.juvM.mu = runif(1,1.4,1.6),  # tight inits here to get it started
                          S.adF.sum.sd = runif(1,0,2), S.adM.sum.sd = runif(1,0,2), 
                          S.adF.win.sd = runif(1,0.2,0.4), S.adM.win.sd = runif(1,0.2,0.4),
                          fecun.juvF.sd = runif(1,0,1), fecun.juvM.sd = runif(1,0,1),   
                          N.adF.spring=c(rpois(1,958.3),rep(NA,19)),
                          N.adM.spring=c(rpois(1,1033.9),rep(NA,19)))}

# Parameters monitored   
parameters <- c("S.adF.sum.mu", "S.adF.win.mu","S.adM.sum.mu", "S.adM.win.mu","fecun.juvF.mu","fecun.juvM.mu", 
                "S.adF.sum.sd", "S.adF.win.sd","S.adM.sum.sd", "S.adM.win.sd","fecun.juvF.sd","fecun.juvM.sd", 
                "N.juvF","N.juvM","N.adM.fall","N.adM.spring", "N.adF.fall","N.adF.spring","S.adF.sum","S.adF.win",
                "S.adM.sum","S.adM.win","fecun.juvF","fecun.juvM",
                "lambda.spring","lambda.fall")

# MCMC settings (use 150,000 ni, 50,000 nb)
# minimal iterations to verify model works
ni <- 150000
nt <- 10
nb <- 50000
nc <- 3

# Call jags from R (bugs model doesn't initiate)
ABDU.Pois_SS.out <- jagsUI(bugs.data, inits, parameters, "ABDU.Pois_SS.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
print(ABDU.Pois_SS.out,digits=3)

########################################################################
# 2-season band recovery model for ABDU
# Female and Male, Brownie
# Juvs transition to adults after first winter, unk start as adult
########################################################################

sink("ABDU.Brownie1.bug")
cat("
    model {
    
    # Priors and constraints, for mean S and f, convert to logit scale
    s.adF.sum.mu ~ dunif(0, 4)          
    s.juvF.win.mu ~ dunif(0, 4)        
    s.adF.win.mu ~ dunif(0, 4)         
    f.juvF.mu ~ dunif(-4, -2)          
    f.adF.mu ~ dunif(-4, -2)        
    s.adM.sum.mu ~ dunif(0, 4)          
    s.juvM.win.mu ~ dunif(0, 4)        
    s.adM.win.mu ~ dunif(0, 4)         
    f.juvM.mu ~ dunif(-4, -2)          
    f.adM.mu ~ dunif(-4, -2)        
    
    # priors for annual SD of survival and reporting rates (logit scale)
    s.adF.sum.sd ~ dunif(0, 1.5)       
    s.juvF.win.sd ~ dunif(0, 1.5)       
    s.adF.win.sd ~ dunif(0, 1.5)        
    f.juvF.sd ~ dunif(0, 0.5)  
    f.adF.sd ~ dunif(0, 0.5)       
    s.adM.sum.sd ~ dunif(0, 1.5)       
    s.juvM.win.sd ~ dunif(0, 1.5)       
    s.adM.win.sd ~ dunif(0, 1.5)        
    f.juvM.sd ~ dunif(0, 0.5)  
    f.adM.sd ~ dunif(0, 0.5)       
    
    # express variances as precision (1/SD^2)
    s.adF.sum.tau <- pow(s.adF.sum.sd,-2)
    s.juvF.win.tau <- pow(s.juvF.win.sd,-2)
    s.adF.win.tau <- pow(s.adF.win.sd,-2)
    f.juvF.tau <- pow(f.juvF.sd,-2)
    f.adF.tau <- pow(f.adF.sd,-2)
    s.adM.sum.tau <- pow(s.adM.sum.sd,-2)
    s.juvM.win.tau <- pow(s.juvM.win.sd,-2)
    s.adM.win.tau <- pow(s.adM.win.sd,-2)
    f.juvM.tau <- pow(f.juvM.sd,-2)
    f.adM.tau <- pow(f.adM.sd,-2)
    
    # Generate seasonal parameters for Survival (juvS is winter only)
    for (y in 1:yrs){
    epsilon.s.juvF.win[y] ~ dnorm(0,s.juvF.win.tau)
    epsilon.s.adF.sum[y] ~ dnorm(0,s.adF.sum.tau)
    epsilon.s.adF.win[y] ~ dnorm(0,s.adF.win.tau)
    epsilon.f.juvF[y] ~ dnorm(0,f.juvF.tau)
    epsilon.f.adF[y] ~ dnorm(0,f.adF.tau)
    epsilon.s.juvM.win[y] ~ dnorm(0,s.juvM.win.tau)
    epsilon.s.adM.sum[y] ~ dnorm(0,s.adM.sum.tau)
    epsilon.s.adM.win[y] ~ dnorm(0,s.adM.win.tau)
    epsilon.f.juvM[y] ~ dnorm(0,f.juvM.tau)
    epsilon.f.adM[y] ~ dnorm(0,f.adM.tau)
    
    logit(s.juvF.win[y]) <- s.juvF.win.mu + epsilon.s.juvF.win[y]
    logit(s.adF.sum[y]) <- s.adF.sum.mu + epsilon.s.adF.sum[y]
    logit(s.adF.win[y]) <- s.adF.win.mu + epsilon.s.adF.win[y]
    logit(f.juvF[y]) <- f.juvF.mu + epsilon.f.juvF[y]
    logit(f.adF[y]) <- f.adF.mu + epsilon.f.adF[y]
    logit(s.juvM.win[y]) <- s.juvM.win.mu + epsilon.s.juvM.win[y]
    logit(s.adM.sum[y]) <- s.adM.sum.mu + epsilon.s.adM.sum[y]
    logit(s.adM.win[y]) <- s.adM.win.mu + epsilon.s.adM.win[y]
    logit(f.juvM[y]) <- f.juvM.mu + epsilon.f.juvM[y]
    logit(f.adM[y]) <- f.adM.mu + epsilon.f.adM[y]
    }
    s.adF.sum[21] <- 1 # dummy values for easier coding
    s.adM.sum[21] <- 1
    
    # Generate multinomial likelihoods for m-arrays
    for (t in 1:(2*yrs)){
    marr.juvF.shot[t,1:(yrs+1)] ~ dmulti(pr.juvF.shot[t,], rel.juvF[t])
    marr.ad_unkF.shot[t,1:(yrs+1)] ~ dmulti(pr.adF.shot[t,], rel.adF[t])
    marr.juvM.shot[t,1:(yrs+1)] ~ dmulti(pr.juvM.shot[t,], rel.juvM[t])
    marr.ad_unkM.shot[t,1:(yrs+1)] ~ dmulti(pr.adM.shot[t,], rel.adM[t])
    }
    
    # Define cell probabilities of the m-arrays
    for (y in 1:yrs){   
    # model prob direct recovery in season [y] (winter banded birds must survive the summer (at adult rates))
    pr.juvF.shot[2*y-1,y] <- s.adF.sum[y]*f.adF[y]       #yearlings recovered at adult rate
    pr.juvF.shot[2*y,y] <- 1*f.juvF[y]                   #fall banded recovered at juvenile rate (fixed to 1)
    pr.adF.shot[2*y-1,y] <- s.adF.sum[y]*f.adF[y]
    pr.adF.shot[2*y,y] <- 1*f.adF[y]
    pr.juvM.shot[2*y-1,y] <- s.adM.sum[y]*f.adM[y]
    pr.juvM.shot[2*y,y] <- 1*f.juvM[y]
    pr.adM.shot[2*y-1,y] <- s.adM.sum[y]*f.adM[y]
    pr.adM.shot[2*y,y] <- 1*f.adM[y]
    
    # probability of surviving to start of second hunting season (2*summer+1winter if winter banded, sum & winter otherwise)
    surv.juvF[2*y-1,y] <- s.adF.sum[y]*s.adF.win[y]*s.adF.sum[y+1]  
    surv.juvF[2*y,y] <- s.juvF.win[y]*s.adF.sum[y+1]
    surv.adF[2*y-1,y] <- s.adF.sum[y]*s.adF.win[y]*s.adF.sum[y+1]
    surv.adF[2*y,y] <- s.adF.win[y]*s.adF.sum[y+1]
    surv.juvM[2*y-1,y] <- s.adM.sum[y]*s.adM.win[y]*s.adM.sum[y+1]
    surv.juvM[2*y,y] <- s.juvM.win[y]*s.adM.sum[y+1]
    surv.adM[2*y-1,y] <- s.adM.sum[y]*s.adM.win[y]*s.adM.sum[y+1]
    surv.adM[2*y,y] <- s.adM.win[y]*s.adM.sum[y+1]
    
    # second and subsequent diagonals
    for (k in (y+1):yrs){
    # model prob recovery = prob survives until k, recovered in k 
    pr.juvF.shot[2*y-1,k] <- surv.juvF[2*y-1,k-1]*f.adF[k]
    pr.juvF.shot[2*y,k] <- surv.juvF[2*y,k-1]*f.adF[k]
    pr.adF.shot[2*y-1,k] <- surv.adF[2*y-1,k-1]*f.adF[k]
    pr.adF.shot[2*y,k] <- surv.adF[2*y,k-1]*f.adF[k]
    pr.juvM.shot[2*y-1,k] <- surv.juvM[2*y-1,k-1]*f.adM[k]
    pr.juvM.shot[2*y,k] <- surv.juvM[2*y,k-1]*f.adM[k]
    pr.adM.shot[2*y-1,k] <- surv.adM[2*y-1,k-1]*f.adM[k]
    pr.adM.shot[2*y,k] <- surv.adM[2*y,k-1]*f.adM[k]
    
    # survival to next hunting season (all adult survival now)
    surv.juvF[2*y-1,k] <- surv.juvF[2*y-1,k-1]*s.adF.win[k-1]*s.adF.sum[k]
    surv.juvF[2*y,k] <- surv.juvF[2*y,k-1]*s.adF.win[k-1]*s.adF.sum[k]
    surv.adF[2*y-1,k] <- surv.adF[2*y-1,k-1]*s.adF.win[k-1]*s.adF.sum[k]
    surv.adF[2*y,k] <- surv.adF[2*y,k-1]*s.adF.win[k-1]*s.adF.sum[k]
    surv.juvM[2*y-1,k] <- surv.juvM[2*y-1,k-1]*s.adM.win[k-1]*s.adM.sum[k]
    surv.juvM[2*y,k] <- surv.juvM[2*y,k-1]*s.adM.win[k-1]*s.adM.sum[k]
    surv.adM[2*y-1,k] <- surv.adM[2*y-1,k-1]*s.adM.win[k-1]*s.adM.sum[k]
    surv.adM[2*y,k] <- surv.adM[2*y,k-1]*s.adM.win[k-1]*s.adM.sum[k]
    } #k
    
    # Left of main diagonal (prob of encounter prior to banding = 0)
    for (l in 1:(y-1)){
    pr.juvF.shot[2*y-1,l] <- 0
    pr.juvF.shot[2*y,l] <- 0
    pr.adF.shot[2*y-1,l] <- 0
    pr.adF.shot[2*y,l] <- 0
    pr.juvM.shot[2*y-1,l] <- 0
    pr.juvM.shot[2*y,l] <- 0
    pr.adM.shot[2*y-1,l] <- 0
    pr.adM.shot[2*y,l] <- 0
    } #l
    } #t
    
    # Last column: probability of non-recovery
    for (t in 1:(2*yrs)){
    pr.juvF.shot[t,(yrs+1)] <- 1-sum(pr.juvF.shot[t,1:yrs])
    pr.adF.shot[t,(yrs+1)] <- 1-sum(pr.adF.shot[t,1:yrs])
    pr.juvM.shot[t,(yrs+1)] <- 1-sum(pr.juvM.shot[t,1:yrs])
    pr.adM.shot[t,(yrs+1)] <- 1-sum(pr.adM.shot[t,1:yrs])
    } #t
    
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data, pool adult and unk m-arrays since not modeling unique parms on uk
bugs.data <- list(yrs=20, marr.juvF.shot=marr.juvF.shot, marr.ad_unkF.shot=marr.ad_unkF.shot,
                  marr.juvM.shot=marr.juvM.shot, marr.ad_unkM.shot=marr.ad_unkM.shot, 
                  rel.juvF = banded[,1], rel.adF = banded[,2]+banded[,3], 
                  rel.juvM = banded[,5], rel.adM = banded[,6]+banded[,7])

# Initial values
inits <- function(){list(s.juvF.win.mu=runif(1,0,2), s.adF.sum.mu=runif(1,0,4), 
                         s.adF.win.mu=runif(1,0,2), f.juvF.mu=runif(1,-4,-2),
                         f.adF.mu=runif(1,-4,-2), s.juvF.win.sd=runif(1,0,1),
                         s.adF.sum.sd=runif(1,0,1), s.adF.win.sd=runif(1,0,1),
                         f.juvF.sd=runif(1,0,0.5), f.adF.sd=runif(1,0,0.5),
                         s.juvM.win.mu=runif(1,0,2), s.adM.sum.mu=runif(1,0,4), 
                         s.adM.win.mu=runif(1,0,2), f.juvM.mu=runif(1,-4,-2),
                         f.adM.mu=runif(1,-4,-2), s.juvM.win.sd=runif(1,0,1),
                         s.adM.sum.sd=runif(1,0,1), s.adM.win.sd=runif(1,0,1),
                         f.juvM.sd=runif(1,0,0.5), f.adM.sd=runif(1,0,0.5))}  

# Parameters monitored
parameters <- c("s.juvF.win.mu", "s.adF.sum.mu", "s.adF.win.mu",
                "f.juvF.mu", "f.adF.mu", "s.juvF.win.sd", "s.adF.sum.sd",
                "s.adF.win.sd", "f.juvF.sd", "f.adF.sd", "s.juvF.win",
                "s.adF.sum", "s.adF.win", "f.juvF", "f.adF",
                "s.juvM.win.mu", "s.adM.sum.mu", "s.adM.win.mu",
                "f.juvM.mu", "f.adM.mu", "s.juvM.win.sd", "s.adM.sum.sd",
                "s.adM.win.sd", "f.juvM.sd", "f.adM.sd", "s.juvM.win",
                "s.adM.sum", "s.adM.win", "f.juvM", "f.adM")

# MCMC settings (31145 sec for 75,000 w 5000 burn on 2 chains, about 9 hours)
## run once more at 150,000 ni, 50,000 nb
# 10,000 iterations, 4902 sec, still hasn't converged on S 
ni <- 150000
nt <- 10
nb <- 50000
nc <- 3

# Call WinBUGS from R #jags version doesn't work, hangs on rel.juvF[1]
#ABDU.Brownie1.bug <- bugs(bugs.data, inits, parameters, "ABDU.Brownie1.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
#print(ABDU.Brownie1.bug, digits = 5)

ABDU.Brownie1.jags <- jagsUI(bugs.data, inits, parameters, "ABDU.Brownie1.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE)
print(ABDU.Brownie1.jags, digits = 5)

########################################################################
# 2-season band recovery model for ABDU
#
# Add bells and whistles, correlations among survival and recovery rates
########################################################################

sink("ABDU.Brownie2.bug")
cat("
    model {
    
    # Priors and constraints
    # Mean survival: 1 JF, 2 JM, 3 AF Sum, 4 AM Sum, 5 AF Winter, 6 AM Winter
    for (s in 1:6){ 
    s.x[s] ~ dunif(0,1)                   # uniform realistic prior on real scale
    s.mu[s] <- logit(s.x[s])              # transform to logit scale
    } # close loop
    
    # Mean Brownie recovery: 1 JF, 2 JM, 3 AF, 4 AM 
    for (f in 1:4){ 
    f.x[f] ~ dunif(0,1)
    f.mu[f] <- logit(f.x[f])
    } # close loop
    
    # Generate correlated variances between males and females for each variable pair
    Omega.js[1:2,1:2] ~ dwish(R.js[,], df.js)  # juvenile survival, 1 F, 2 M
    Sigma.js[1:2,1:2] <- inverse(Omega.js[,])
    Omega.ss[1:2,1:2] ~ dwish(R.ss[,], df.ss)  # adult summer survival
    Sigma.ss[1:2,1:2] <- inverse(Omega.ss[,])
    Omega.ws[1:2,1:2] ~ dwish(R.ws[,], df.ws)  # adult winter survival
    Sigma.ws[1:2,1:2] <- inverse(Omega.ws[,])
    Omega.jf[1:2,1:2] ~ dwish(R.jf[,], df.jf) # juvenile recovery 
    Sigma.jf[1:2,1:2] <- inverse(Omega.jf[,])
    Omega.af[1:2,1:2] ~ dwish(R.af[,], df.af) # adult recovery rate
    Sigma.af[1:2,1:2] <- inverse(Omega.af[,])
    
    # Generate seasonal parameters for Survival and recovery (1 is F, 2 is M)
    for (y in 1:yrs){
    eps.js[1:2,y] ~ dmnorm(zero.js[], Omega.js[,])
    eps.ss[1:2,y] ~ dmnorm(zero.ss[], Omega.ss[,])
    eps.ws[1:2,y] ~ dmnorm(zero.ws[], Omega.ws[,])
    eps.jf[1:2,y] ~ dmnorm(zero.jf[], Omega.jf[,])
    eps.af[1:2,y] ~ dmnorm(zero.af[], Omega.af[,])
    for (s in 1:2){
    logit(js[s,y]) <- s.mu[s] + eps.js[s,y]
    logit(ss[s,y]) <- s.mu[s+2] + eps.ss[s,y]
    logit(ws[s,y]) <- s.mu[s+4] + eps.ws[s,y]
    logit(jf[s,y]) <- f.mu[s] + eps.jf[s,y]
    logit(af[s,y]) <- f.mu[s+2] + eps.af[s,y]}
    } # close y
    ss[1,21] <- 1 # dummy values for extra summer survival for easier coding
    ss[2,21] <- 1
    
    # Generate multinomial likelihoods for m-arrays
    for (t in 1:(2*yrs)){
    marr.juvF.shot[t,1:(yrs+1)] ~ dmulti(pr.juvF.shot[t,], rel.juvF[t])
    marr.ad_unkF.shot[t,1:(yrs+1)] ~ dmulti(pr.adF.shot[t,], rel.adF[t])
    marr.juvM.shot[t,1:(yrs+1)] ~ dmulti(pr.juvM.shot[t,], rel.juvM[t])
    marr.ad_unkM.shot[t,1:(yrs+1)] ~ dmulti(pr.adM.shot[t,], rel.adM[t])
    }
    
    # Define cell probabilities of the m-arrays
    for (y in 1:yrs){   
    # model prob direct recovery in season [y] 
    pr.juvF.shot[2*y-1,y] <- ss[1,y] * af[1,y]   # juv marked in spring survives summer as adult, recovered as adult
    pr.juvF.shot[2*y,y] <- 1 * jf[1,y]
    pr.adF.shot[2*y-1,y] <- ss[1,y] * af[1,y]
    pr.adF.shot[2*y,y] <- 1 * af[1,y]
    pr.juvM.shot[2*y-1,y] <- ss[2,y] * af[2,y]
    pr.juvM.shot[2*y,y] <- 1 * jf[2,y]
    pr.adM.shot[2*y-1,y] <- ss[2,y] * af[2,y]
    pr.adM.shot[2*y,y] <- 1 * af[2,y]
    
    # probability of surviving to start of second hunting season (2*summer+1winter if winter banded, sum & winter otherwise)
    surv.jF[2*y-1,y] <- ss[1,y] * ws[1,y] * ss[1,y+1]  
    surv.jF[2*y,y] <- js[1,y] * ss[1,y+1]
    surv.aF[2*y-1,y] <- ss[1,y] * ws[1,y] * ss[1,y+1]
    surv.aF[2*y,y] <- ws[1,y] * ss[1,y+1]
    surv.jM[2*y-1,y] <- ss[2,y] * ws[2,y] * ss[2,y+1]
    surv.jM[2*y,y] <- js[2,y] * ss[2,y+1]
    surv.aM[2*y-1,y] <- ss[2,y] * ws[2,y] * ss[2,y+1]
    surv.aM[2*y,y] <- ws[2,y] * ss[2,y+1]
    
    # second and subsequent diagonals
    for (k in (y+1):yrs){
    # model prob recovery = prob survives until k, recovered in k 
    pr.juvF.shot[2*y-1,k] <- surv.jF[2*y-1,k-1] * af[1,k]
    pr.juvF.shot[2*y,k] <- surv.jF[2*y,k-1] * af[1,k]
    pr.adF.shot[2*y-1,k] <- surv.aF[2*y-1,k-1] * af[1,k]
    pr.adF.shot[2*y,k] <- surv.aF[2*y,k-1] * af[1,k]
    pr.juvM.shot[2*y-1,k] <- surv.jM[2*y-1,k-1] * af[2,k]
    pr.juvM.shot[2*y,k] <- surv.jM[2*y,k-1] * af[2,k]
    pr.adM.shot[2*y-1,k] <- surv.aM[2*y-1,k-1] * af[2,k]
    pr.adM.shot[2*y,k] <- surv.aM[2*y,k-1] * af[2,k]
    
    # survival to next hunting season (all adult survival now)
    surv.jF[2*y-1,k] <- surv.jF[2*y-1,k-1] * ss[1,k-1] * ws[1,k]
    surv.jF[2*y,k] <- surv.jF[2*y,k-1] * ss[1,k-1] * ws[1,k]
    surv.aF[2*y-1,k] <- surv.aF[2*y-1,k-1] * ss[1,k-1] * ws[1,k]
    surv.aF[2*y,k] <- surv.aF[2*y,k-1] * ss[1,k-1] * ws[1,k]
    surv.jM[2*y-1,k] <- surv.jM[2*y-1,k-1] * ss[2,k-1] * ws[2,k]
    surv.jM[2*y,k] <- surv.jM[2*y,k-1] * ss[2,k-1] * ws[2,k]
    surv.aM[2*y-1,k] <- surv.aM[2*y-1,k-1] * ss[2,k-1] * ws[2,k]
    surv.aM[2*y,k] <- surv.aM[2*y,k-1] * ss[2,k-1] * ws[2,k]
    } #k
    
    # Left of main diagonal (prob of encounter prior to banding = 0)
    for (l in 1:(y-1)){
    pr.juvF.shot[2*y-1,l] <- 0
    pr.juvF.shot[2*y,l] <- 0
    pr.adF.shot[2*y-1,l] <- 0
    pr.adF.shot[2*y,l] <- 0
    pr.juvM.shot[2*y-1,l] <- 0
    pr.juvM.shot[2*y,l] <- 0
    pr.adM.shot[2*y-1,l] <- 0
    pr.adM.shot[2*y,l] <- 0
    } #l
    } #y
    
    # Last column: probability of non-recovery
    for (t in 1:(2*yrs)){
    pr.juvF.shot[t,(yrs+1)] <- 1-sum(pr.juvF.shot[t,1:yrs])
    pr.adF.shot[t,(yrs+1)] <- 1-sum(pr.adF.shot[t,1:yrs])
    pr.juvM.shot[t,(yrs+1)] <- 1-sum(pr.juvM.shot[t,1:yrs])
    pr.adM.shot[t,(yrs+1)] <- 1-sum(pr.adM.shot[t,1:yrs])
    } #t
    
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data, pool adult and unk m-arrays since not modeling unique parms on uk
bugs.data <- list(yrs=20, marr.juvF.shot=marr.juvF.shot, marr.ad_unkF.shot=marr.ad_unkF.shot,
                  marr.juvM.shot=marr.juvM.shot, marr.ad_unkM.shot=marr.ad_unkM.shot, 
                  rel.juvF = banded[,1], rel.adF = banded[,2]+banded[,3], 
                  rel.juvM = banded[,5], rel.adM = banded[,6]+banded[,7],
                  zero.js = rep(0,2), df.js = 18, R.js = diag(2),
                  zero.ss = rep(0,2), df.ss = 18, R.ss = diag(2),
                  zero.ws = rep(0,2), df.ws = 18, R.ws = diag(2),
                  zero.jf = rep(0,2), df.jf = 18, R.jf = diag(2),
                  zero.af = rep(0,2), df.af = 18, R.af = diag(2))

# Initial values
inits <- function(){list(s.x=c(rep(runif(1,0,1),6)), f.x=c(rep(runif(1,0,1),4)))}  

# Parameters monitored
parameters <- c("s.x", "f.x", "Sigma.js", "Sigma.ss",  "Sigma.ws", "Sigma.jf", "Sigma.af", 
                "js", "ss","ws", "jf", "af")

# MCMC settings (150000, 50000nb achieves convergence on all parameters)
ni <- 15000
nt <- 10
nb <- 5000
nc <- 1

# Call jags from R
# can hang at initiation (with delayed error message), monitor for first few minutes if running in parallel
ABDU.Brownie2.jags <- jagsUI(bugs.data, inits, parameters, "ABDU.Brownie2.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=FALSE)
print(ABDU.Brownie2.jags, digits = 4)

#######################
# Part 3: Generate an IPM using Lincoln estimates for N and F, band recoveries for S
########################

sink("ABDU.IPM.bug")
cat("
    model {
    
    #---------------------------------------------
    # 1. Priors and constraints
    #---------------------------------------------
    
    # observed initial Lincoln estimates/1000 (N) for combined juv + adult females (af) and males (am)
    N.f.Feb[1] ~ dpois(958.3)    
    N.m.Feb[1] ~ dpois(1033.9)
    
    # uniform 0-1 priors for logit survival (s) and Brownie recovery (f) rates (Koons & Schaub 2015)
    s.af.sum.mu ~ dunif(0,1)                            
    ls.af.sum.mu <- log(s.af.sum.mu/(1-s.af.sum.mu))    
    s.am.sum.mu ~ dunif(0,1)                            
    ls.am.sum.mu <- log(s.am.sum.mu/(1-s.am.sum.mu))    
    s.af.win.mu ~ dunif(0,1)                            
    ls.af.win.mu <- log(s.af.win.mu/(1-s.af.win.mu))    
    s.am.win.mu ~ dunif(0,1)                            
    ls.am.win.mu <- log(s.am.win.mu/(1-s.am.win.mu))    
    s.jf.sum.mu ~ dunif(0,1)                            # Give yrlg separate S in 1st summer?
    ls.jf.sum.mu <- log(s.jf.sum.mu/(1-s.jf.sum.mu))    
    s.jm.sum.mu ~ dunif(0,1)                            
    ls.jm.sum.mu <- log(s.jm.sum.mu/(1-s.jm.sum.mu))    
    s.jf.win.mu ~ dunif(0,1)                            
    ls.jf.win.mu <- log(s.jf.win.mu/(1-s.jf.win.mu))    
    s.jm.win.mu ~ dunif(0,1)                            
    ls.jm.win.mu <- log(s.jm.win.mu/(1-s.jm.win.mu))    
    
    f.af.mu ~ dunif(0,0.3)                            
    lf.af.mu <- log(f.af.mu/(1-f.af.mu))    
    f.am.mu ~ dunif(0,0.3)                            
    lf.am.mu <- log(f.am.mu/(1-f.am.mu))    
    f.jf.mu ~ dunif(0,0.3)                            
    lf.jf.mu <- log(f.jf.mu/(1-f.jf.mu))    
    f.jm.mu ~ dunif(0,0.3)                            
    lf.jm.mu <- log(f.jm.mu/(1-f.jm.mu))    
    
    #### Revisit this, F with added parameter for sex ratio at fledge
    F.jf.mu ~ dunif(0.5,2.5)      # prior for mean fecundity rates (F sex spec because fall sex ratio unbalanced)
    F.jm.mu ~ dunif(0.5,2.5)      
    
    # annual sd, logit scale (from preliminary analyses in MARK, these are all < 0.5)
    s.af.sum.sd ~ dunif(0.05,2)                            
    s.am.sum.sd ~ dunif(0.05,2)                            
    s.af.win.sd ~ dunif(0.05,2)                            
    s.am.win.sd ~ dunif(0.05,2)                            
    s.jf.sum.sd ~ dunif(0.05,2)                            
    s.jm.sum.sd ~ dunif(0.05,2)                            
    s.jf.win.sd ~ dunif(0.05,2)                            
    s.jm.win.sd ~ dunif(0.05,2)                            
    
    f.af.sd ~ dunif(0.05,2)                            
    f.am.sd ~ dunif(0.05,2)                            
    f.jf.sd ~ dunif(0.05,2)                            
    f.jm.sd ~ dunif(0.05,2)                            
    
    F.jf.sd ~ dunif(0.05,2)    
    F.jm.sd ~ dunif(0.05,2)      
    
    # annual sd expressed as precision, tau = (1/sd^2)
    s.af.sum.tau <- pow(s.af.sum.sd,-2)                            
    s.am.sum.tau <- pow(s.am.sum.sd,-2)                            
    s.af.win.tau <- pow(s.af.win.sd,-2)                            
    s.am.win.tau <- pow(s.am.win.sd,-2)                            
    s.jf.sum.tau <- pow(s.jf.sum.sd,-2)                            
    s.jm.sum.tau <- pow(s.jm.sum.sd,-2)                            
    s.jf.win.tau <- pow(s.jf.win.sd,-2)                            
    s.jm.win.tau <- pow(s.jm.win.sd,-2)                            
    
    f.af.tau <- pow(f.af.sd,-2)                            
    f.am.tau <- pow(f.am.sd,-2)                            
    f.jf.tau <- pow(f.jf.sd,-2)                            
    f.jm.tau <- pow(f.jm.sd,-2)                            
    
    F.jf.tau <- pow(F.jf.sd,-2)    
    F.jm.tau <- pow(F.jm.sd,-2)     
    
    #---------------------------------------------
    # 2. Generate annual parameter estimates (S, f, F)
    #---------------------------------------------
    
    for (y in 1:yrs){
    logit.s.af.sum[y] ~ dnorm(ls.af.sum.mu,s.af.sum.tau)
    logit(s.af.sum[y]) <- logit.s.af.sum[y] 
    logit.s.am.sum[y] ~ dnorm(ls.am.sum.mu,s.am.sum.tau)
    logit(s.am.sum[y]) <- logit.s.am.sum[y] 
    logit.s.af.win[y] ~ dnorm(ls.af.win.mu,s.af.win.tau)
    logit(s.af.win[y]) <- logit.s.af.win[y] 
    logit.s.am.win[y] ~ dnorm(ls.am.win.mu,s.am.win.tau)
    logit(s.am.win[y]) <- logit.s.am.win[y] 
    logit.s.jf.sum[y] ~ dnorm(ls.jf.sum.mu,s.jf.sum.tau)
    logit(s.jf.sum[y]) <- logit.s.jf.sum[y] 
    logit.s.jm.sum[y] ~ dnorm(ls.jm.sum.mu,s.jm.sum.tau)
    logit(s.jm.sum[y]) <- logit.s.jm.sum[y] 
    logit.s.jf.win[y] ~ dnorm(ls.jf.win.mu,s.jf.win.tau)
    logit(s.jf.win[y]) <- logit.s.jf.win[y] 
    logit.s.jm.win[y] ~ dnorm(ls.jm.win.mu,s.jm.win.tau)
    logit(s.jm.win[y]) <- logit.s.jm.win[y] 
    logit.f.af[y] ~ dnorm(lf.af.mu,f.af.tau)
    logit(f.af[y]) <- logit.f.af[y] 
    logit.f.am[y] ~ dnorm(lf.am.mu,f.am.tau)
    logit(f.am[y]) <- logit.f.am[y] 
    logit.f.jf[y] ~ dnorm(lf.jf.mu,f.jf.tau)
    logit(f.jf[y]) <- logit.f.jf[y] 
    logit.f.jm[y] ~ dnorm(lf.jm.mu,f.jm.tau)
    logit(f.jm[y]) <- logit.f.jm[y] 
    Fecun.jf[y] ~ dnorm(F.jf.mu,F.jf.tau) 
    F.jf[y] <-max(0,Fecun.jf[y])                # avoiding I(0,) format
    Fecun.jm[y] ~ dnorm(F.jm.mu,F.jm.tau) 
    F.jm[y] <-max(0,Fecun.jm[y])
    }
    s.af.sum[21] <- 1 # dummy variables for parameter references exceeding matrix boundaries
    s.am.sum[21] <- 1
    s.jf.sum[21] <- 1
    s.jm.sum[21] <- 1
    
    #-------------------------------------------
    # 3a. Likelihood for population components from Lincoln estimates
    #-------------------------------------------
    # millions of birds, so ignoring demographic stochasticity
    # apportion spring 1 population estimates according to long-term projected age ratios
    N.af.Feb[1] <- N.f.Feb[1] * (1-0.57)  
    N.am.Feb[1] <- N.m.Feb[1] * (1-0.42)  
    N.jf.Feb[1] <- N.f.Feb[1] * 0.57
    N.jm.Feb[1] <- N.f.Feb[1] * 0.42
    
    
    # estimate fall population sizes in year 1
    N.af.Aug[1] <- N.af.Feb[1] * s.af.sum[1]  
    N.am.Aug[1] <- N.am.Feb[1] * s.am.sum[1]  
    N.jf.Aug[1] <- N.af.Aug[1] * F.jf[1]
    N.jm.Aug[1] <- N.af.Aug[1] * F.jm[1]
    
    for (t in 2:yrs){
    N.af.Feb[t] <- N.af.Aug[t-1] * s.af.win[t-1]
    N.am.Feb[t] <- N.am.Aug[t-1] * s.am.win[t-1]
    N.jf.Feb[t] <- N.jf.Aug[t-1] * s.jf.win[t-1]
    N.jm.Feb[t] <- N.jm.Aug[t-1] * s.jm.win[t-1]
    # combined ad + juv population is what Lincoln (and BPOP) estimates
    N.f.Feb[t] <- N.af.Feb[t] + N.jf.Feb[t]        
    N.m.Feb[t] <- N.am.Feb[t] + N.jm.Feb[t]        
    # juveniles (yearlings) graduate to adults in this step
    N.af.Aug[t] <- N.af.Feb[t] * s.af.sum[t]  + N.jf.Feb[t] * s.jf.sum[t]
    N.am.Aug[t] <- N.am.Feb[t] * s.am.sum[t]  + N.jm.Feb[t] * s.jm.sum[t]
    N.jf.Aug[t] <- N.af.Aug[t] * F.jf[t]
    N.jm.Aug[t] <- N.af.Aug[t] * F.jm[t]
    }
    
    # Derived parms, more efficient to estimate from saved posteriors unless I want intrinsic correlations
    #  lambda.BPOP[t-1] <- (N.f.Feb[t]+N.m.Feb[t])/(N.f.Feb[t-1]+N.m.Feb[t-1])
    #  lambda.Fall[t-1] <- (N.af.Aug[t]+N.am.Aug[t]+N.jf.Aug[t]+N.jm.Aug[t])/(N.af.Aug[t-1]+N.am.Aug[t-1]+N.jf.Aug[t-1]+N.jm.Aug[t-1])
    #  sexratio.BPOP[t] <- (N.m.Feb[t])/(N.f.Feb[t]+N.m.Feb[t])
    
    # Observation process, population counts
    for (t in 1:yrs) {
    jf.obs.Aug[t] ~ dpois(N.jf.Aug[t]) 
    jm.obs.Aug[t] ~ dpois(N.jm.Aug[t]) 
    af.obs.Aug[t] ~ dpois(N.af.Aug[t]) 
    am.obs.Aug[t] ~ dpois(N.am.Aug[t]) 
    f.obs.Feb[t] ~ dpois(N.f.Feb[t]) 
    m.obs.Feb[t] ~ dpois(N.m.Feb[t]) 
    }
    
    #-------------------------------------------
    # 3a. Likelihood for band recovery data (S and f)
    #-------------------------------------------
    
    # Generate multinomial likelihoods for m-arrays (modified from Kery & Schaub)
    # 2 rows for each year, Feb bandings, then Aug bandings
    # also have m-arrays for other dead encounters and live recaps, but they add little precision, shot only here
    
    for (t in 1:(2*yrs)){
    marr.juvF.shot[t,1:(yrs+1)] ~ dmulti(pr.jf.shot[t,], rel.jf[t])
    marr.ad_unkF.shot[t,1:(yrs+1)] ~ dmulti(pr.af.shot[t,], rel.af[t])
    marr.juvM.shot[t,1:(yrs+1)] ~ dmulti(pr.jm.shot[t,], rel.jm[t])
    marr.ad_unkM.shot[t,1:(yrs+1)] ~ dmulti(pr.am.shot[t,], rel.am[t])
    }
    
    # Define cell probabilities of the m-arrays
    for (y in 1:yrs){   
    # model prob direct recovery (postseason juvs survive sum as juvs, but harvested as adults)
    pr.jf.shot[2*y-1,y] <- s.jf.sum[y] * f.af[y]
    pr.jf.shot[2*y,y] <- 1 * f.jf[y]              # 1  denotes all fall-banded birds survived summer
    pr.af.shot[2*y-1,y] <- s.af.sum[y] * f.af[y]
    pr.af.shot[2*y,y] <- 1 * f.af[y]
    pr.jm.shot[2*y-1,y] <- s.jm.sum[y] * f.am[y]
    pr.jm.shot[2*y,y] <- 1 * f.jm[y]
    pr.am.shot[2*y-1,y] <- s.am.sum[y] * f.am[y]
    pr.am.shot[2*y,y] <- 1 * f.am[y]
    
    # cumulative probability of surviving to start of 2nd hunting season (surv), seasonal components (s)
    surv.jf[2*y-1,y] <- s.jf.sum[y] * s.af.win[y] * s.af.sum[y+1] 
    surv.jf[2*y,y] <- s.jf.win[y] * s.jf.sum[y+1]
    surv.af[2*y-1,y] <- s.af.sum[y] * s.af.win[y] * s.af.sum[y+1]
    surv.af[2*y,y] <- s.af.win[y] * s.af.sum[y+1]
    surv.jm[2*y-1,y] <- s.jm.sum[y] * s.am.win[y] * s.am.sum[y+1]
    surv.jm[2*y,y] <- s.jm.win[y] * s.jm.sum[y+1]
    surv.am[2*y-1,y] <- s.am.sum[y] * s.am.win[y] * s.am.sum[y+1]
    surv.am[2*y,y] <- s.am.win[y] * s.am.sum[y+1]
    } # y for first diagonal
    
    # second y loop 
    for (y in 1:yrs){   
    for (j in (y+1):yrs){
    # model prob recovery = prob survives until j, recovered in j 
    pr.jf.shot[2*y-1,j] <- surv.jf[2*y-1,j-1]*f.af[j]
    pr.jf.shot[2*y,j] <- surv.jf[2*y,j-1]*f.af[j]
    pr.af.shot[2*y-1,j] <- surv.af[2*y-1,j-1]*f.af[j]
    pr.af.shot[2*y,j] <- surv.af[2*y,j-1]*f.af[j]
    pr.jm.shot[2*y-1,j] <- surv.jm[2*y-1,j-1]*f.am[j]
    pr.jm.shot[2*y,j] <- surv.jm[2*y,j-1]*f.am[j]
    pr.am.shot[2*y-1,j] <- surv.am[2*y-1,j-1]*f.am[j]
    pr.am.shot[2*y,j] <- surv.am[2*y,j-1]*f.am[j]
    # model survival to next hunting season
    surv.jf[2*y-1,j] <- surv.jf[2*y-1,j-1]*s.af.win[j-1]*s.af.sum[j] 
    surv.jf[2*y,j]   <- surv.jf[2*y,j-1]*s.af.win[j-1]*s.af.sum[j]
    surv.af[2*y-1,j] <- surv.af[2*y-1,j-1]*s.af.win[j-1]*s.af.sum[j]
    surv.af[2*y,j]   <- surv.af[2*y,j-1]*s.af.win[j-1]*s.af.sum[j]
    surv.jm[2*y-1,j] <- surv.jm[2*y-1,j-1]*s.am.win[j-1]*s.am.sum[j] 
    surv.jm[2*y,j]   <- surv.jm[2*y,j-1]*s.am.win[j-1]*s.am.sum[j]
    surv.am[2*y-1,j] <- surv.am[2*y-1,j-1]*s.am.win[j-1]*s.am.sum[j]
    surv.am[2*y,j]   <- surv.am[2*y,j-1]*s.am.win[j-1]*s.am.sum[j]
    } #j
    
    } # end second y loop
    
    # third y loop fills in remainder of dmulti matrix
    for (y in 1:yrs){  
    # Left of main diagonal (prob of encounter prior to banding = 0)
    for (l in 1:(y-1)){
    pr.jf.shot[2*y-1,l] <- 0
    pr.jf.shot[2*y,l] <- 0
    pr.af.shot[2*y-1,l] <- 0
    pr.af.shot[2*y,l] <- 0
    pr.jm.shot[2*y-1,l] <- 0
    pr.jm.shot[2*y,l] <- 0
    pr.am.shot[2*y-1,l] <- 0
    pr.am.shot[2*y,l] <- 0
    } #l
    } #y
    
    # Last column: probability of non-recovery
    for (t in 1:(2*yrs)){
    pr.jf.shot[t,(yrs+1)] <- 1-sum(pr.jf.shot[t,1:yrs])
    pr.af.shot[t,(yrs+1)] <- 1-sum(pr.af.shot[t,1:yrs])
    pr.jm.shot[t,(yrs+1)] <- 1-sum(pr.jm.shot[t,1:yrs])
    pr.am.shot[t,(yrs+1)] <- 1-sum(pr.am.shot[t,1:yrs])
    } #t
    }  # Hail Mary, end bugs model  
    ",fill=TRUE)
sink()

# Bundle data 
bugs.data <- list(yrs=20, jf.obs.Aug=trunc(Lincoln.out[,1]/1000), jm.obs.Aug=trunc(Lincoln.out[,3]/1000), af.obs.Aug=trunc(Lincoln.out[,5]/1000),
                  am.obs.Aug=trunc(Lincoln.out[,7]/1000), f.obs.Feb=trunc(Lincoln.out[,9]/1000), m.obs.Feb=trunc(Lincoln.out[,11]/1000), 
                  marr.juvF.shot=marr.juvF.shot, marr.ad_unkF.shot=marr.ad_unkF.shot, marr.juvM.shot=marr.juvM.shot, marr.ad_unkM.shot=marr.ad_unkM.shot,
                  rel.jf=rowSums(marr.juvF.shot), rel.af=rowSums(marr.ad_unkF.shot),
                  rel.jm=rowSums(marr.juvM.shot), rel.am=rowSums(marr.ad_unkM.shot))

# Initial values (starting chains in target range)  
inits <- function() {list(s.af.sum.mu = runif(1,0,3), s.am.sum.mu = runif(1,0,3), s.jf.sum.mu = runif(1,0,3), s.jm.sum.mu = runif(1,0,3),
                          s.af.win.mu = runif(1,0,3), s.am.win.mu = runif(1,0,3), s.jf.win.mu = runif(1,0,3), s.jm.win.mu = runif(1,0,3),
                          f.af.mu = runif(1,-4,-1), f.am.mu = runif(1,-4,-1), f.jf.mu = runif(1,-4,-1), f.jm.mu = runif(1,-4,-1),
                          F.jf.mu = runif(1,0.5,2), F.jm.mu = runif(1,0.5,2),   
                          N.f.Feb=c(rpois(1,958.3),rep(NA,19)),
                          N.m.Feb=c(rpois(1,1033.9),rep(NA,19)))}

# Parameters monitored   
parameters <- c("N.af.Feb","N.am.Feb","N.jf.Feb","N.jm.Feb","N.f.Feb","N.m.Feb","N.af.Aug","N.am.Aug","N.jf.Aug","N.jm.Aug",
                "s.af.sum.mu","s.am.sum.mu","s.af.win.mu","s.am.win.mu","s.jf.sum.mu","s.jm.sum.mu","s.jf.win.mu","s.jm.win.mu",
                "s.af.sum.sd","s.am.sum.sd","s.af.win.sd","s.am.win.sd","s.jf.sum.sd","s.jm.sum.sd","s.jf.win.sd","s.jm.win.sd",
                "f.af.mu","f.am.mu","f.jf.mu","f.jm.mu","f.af.sd","f.am.sd","f.jf.sd","f.jm.sd",
                "s.af.sum","s.af.win","s.am.sum","s.am.win","s.jf.sum","s.jf.win","s.jm.sum","s.jm.win",
                "f.af","f.am","f.jf","f.jm")

# MCMC settings (just trying to get this thing to go)
ni <- 100
nt <- 1
nb <- 0
nc <- 3

# Call WinBUGS from R (bugs model gives undefined real result, jags hangs on rel.jf[1])
#ABDU.IPM.out <- bugs(bugs.data, inits, parameters, "ABDU.IPM.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,debug=TRUE,bugs.directory=bugs.dir,working.directory=getwd())
#print(ABDU.IPM.out, digits=4)

ABDU.IPM.out <- jagsUI(bugs.data, inits, parameters, "ABDU.IPM.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
#print(ABDU.IPM.out, digits=4)


##################################################################
# simpler version, Females only, juvs transition to adults in Feb
##################################################################

sink("ABDU.IPM2.bug")
cat("
    model {
    
    #---------------------------------------------
    # 1. Priors and constraints
    #---------------------------------------------
    
    # observed initial Lincoln estimates/1000 (N) for combined juv + adult females (af) and males (am)
    N.f.Feb[1] ~ dpois(958.3)    
    
    # uniform 0-1 priors for logit survival (s) and Brownie recovery (f) rates (Koons & Schaub 2015)
    s.af.sum.mu ~ dunif(0,1)                            
    ls.af.sum.mu <- log(s.af.sum.mu/(1-s.af.sum.mu))    
    s.af.win.mu ~ dunif(0,1)                            
    ls.af.win.mu <- log(s.af.win.mu/(1-s.af.win.mu))    
    s.jf.win.mu ~ dunif(0,1)                            
    ls.jf.win.mu <- log(s.jf.win.mu/(1-s.jf.win.mu))    
    
    f.af.mu ~ dunif(0,1)                            
    lf.af.mu <- log(f.af.mu/(1-f.af.mu))    
    f.jf.mu ~ dunif(0,1)                            
    lf.jf.mu <- log(f.jf.mu/(1-f.jf.mu))    
    
    F.jf.mu ~ dunif(0.5,2.5)      # prior for mean fecundity rates (F sex spec because fall sex ratio unbalanced)
    
    # annual sd, logit scale (from preliminary analyses in MARK, these are all < 0.5)
    s.af.sum.sd ~ dunif(0,2)                            
    s.af.win.sd ~ dunif(0,2)                            
    s.jf.win.sd ~ dunif(0,2)                            
    f.af.sd ~ dunif(0,2)                            
    f.jf.sd ~ dunif(0,2)                            
    F.jf.sd ~ dunif(0,2)    
    
    # annual sd expressed as precision, tau = (1/sd^2)
    s.af.sum.tau <- pow(s.af.sum.sd,-2)                            
    s.af.win.tau <- pow(s.af.win.sd,-2)                            
    s.jf.win.tau <- pow(s.jf.win.sd,-2)                            
    f.af.tau <- pow(f.af.sd,-2)                            
    f.jf.tau <- pow(f.jf.sd,-2)                            
    F.jf.tau <- pow(F.jf.sd,-2)    
    
    #---------------------------------------------
    # 2. Generate annual parameter estimates (S, f, F)
    #---------------------------------------------
    
    for (y in 1:yrs){
    logit.s.af.sum[y] ~ dnorm(ls.af.sum.mu,s.af.sum.tau)
    logit(s.af.sum[y]) <- logit.s.af.sum[y] 
    logit.s.af.win[y] ~ dnorm(ls.af.win.mu,s.af.win.tau)
    logit(s.af.win[y]) <- logit.s.af.win[y] 
    logit.s.jf.win[y] ~ dnorm(ls.jf.win.mu,s.jf.win.tau)
    logit(s.jf.win[y]) <- logit.s.jf.win[y] 
    logit.f.af[y] ~ dnorm(lf.af.mu,f.af.tau)
    logit(f.af[y]) <- logit.f.af[y] 
    logit.f.jf[y] ~ dnorm(lf.jf.mu,f.jf.tau)
    logit(f.jf[y]) <- logit.f.jf[y] 
    Fecun.jf[y] ~ dnorm(F.jf.mu,F.jf.tau) 
    F.jf[y] <-max(0,Fecun.jf[y])                # avoiding I(0,) format
    }
    
    #-------------------------------------------
    # 3a. Likelihood for population components from Lincoln estimates
    #-------------------------------------------
    # millions of birds, so ignoring demographic stochasticity
    # estimate fall population sizes in year 1
    N.af.Aug[1] <- N.f.Feb[1] * s.af.sum[1]  
    N.jf.Aug[1] <- N.af.Aug[1] * F.jf[1]
    
    for (t in 2:yrs){
    N.f.Feb[t] <- N.af.Aug[t-1] * s.af.win[t-1] + N.jf.Aug[t-1] * s.jf.win[t-1]
    N.af.Aug[t] <- N.f.Feb[t] * s.af.sum[t]  
    N.jf.Aug[t] <- N.af.Aug[t] * F.jf[t]
    }
    
    # Observation process, population counts
    for (t in 1:yrs) {
    jf.obs.Aug[t] ~ dpois(N.jf.Aug[t]) 
    af.obs.Aug[t] ~ dpois(N.af.Aug[t]) 
    f.obs.Feb[t] ~ dpois(N.f.Feb[t]) 
    }
    
    #-------------------------------------------
    # 3a. Likelihood for band recovery data (S and f)
    #-------------------------------------------
    
    # Generate multinomial likelihoods for m-arrays (modified from Kery & Schaub)
    # 2 rows for each year, Feb bandings, then Aug bandings
    # also have m-arrays for other dead encounters and live recaps, but they add little precision, shot only here
    
    for (t in 1:(2*yrs)){
    marr.juvF.shot[t,1:(yrs+1)] ~ dmulti(pr.jf.shot[t,], rel.jf[t])
    marr.ad_unkF.shot[t,1:(yrs+1)] ~ dmulti(pr.af.shot[t,], rel.af[t])
    }
    
    # Calculate the number of birds released each banding season (releases in all 40 rows)
    for (t in 1:(2*yrs)){
    rel.jf[t] <- sum(marr.juvF.shot[t,])
    rel.af[t] <- sum(marr.ad_unkF.shot[t,])
    }
    
    # Define cell probabilities of the m-arrays
    for (y in 1:yrs){   
    # cumulative probability of surviving to start of hunting season (surv), seasonal components (s)
    surv.jf[2*y-1,y] <- s.af.sum[y]  
    surv.jf[2*y,y] <- 1
    surv.af[2*y-1,y] <- s.af.sum[y]
    surv.af[2*y,y] <- 1
    
    # model prob direct recovery (postseason juvs survive sum as juvs, but harvested as adults)
    pr.jf.shot[2*y-1,y] <- surv.jf[2*y-1,y]*f.af[y]
    pr.jf.shot[2*y,y] <- surv.jf[2*y,y]*f.jf[y]
    pr.af.shot[2*y-1,y] <- surv.af[2*y-1,y]*f.af[y]
    pr.af.shot[2*y,y] <- surv.af[2*y,y]*f.af[y]
    } # y for first diagonal
    
    # second y loop omits last 2 rows and facilitates 2nd diagonal equations
    for (y in 1:(yrs-1)){   
    # second diagonal (wherein preseason juvs experience juvenile survival in winter, postseason juvs are now adults)
    for (j in (y+1):(y+1)){
    surv.jf[2*y-1,j] <- surv.jf[2*y-1,j-1]*s.af.win[j-1]*s.af.sum[j] 
    surv.jf[2*y,j]   <- surv.jf[2*y,j-1]*s.jf.win[j-1]*s.af.sum[j]
    surv.af[2*y-1,j] <- surv.af[2*y-1,j-1]*s.af.win[j-1]*s.af.sum[j]
    surv.af[2*y,j]   <- surv.af[2*y,j-1]*s.af.win[j-1]*s.af.sum[j]
    
    # model prob recovery = prob survives until j, recovered in j 
    pr.jf.shot[2*y-1,j] <- surv.jf[2*y-1,j]*f.af[j]
    pr.jf.shot[2*y,j] <- surv.jf[2*y,j]*f.af[j]
    pr.af.shot[2*y-1,j] <- surv.jf[2*y-1,j]*f.af[j]
    pr.af.shot[2*y,j] <- surv.jf[2*y,j]*f.af[j]
    } #k
    
    # third and subsequent diagonals (all cohorts have transitioned to adulthood)
    for (k in (y+2):yrs){
    surv.jf[2*y-1,k] <- surv.jf[2*y-1,k-1]*s.af.win[k-1]*s.af.sum[k]
    surv.jf[2*y,k] <- surv.jf[2*y,k-1]*s.af.win[k-1]*s.af.sum[k]
    surv.af[2*y-1,k] <- surv.af[2*y-1,k-1]*s.af.win[k-1]*s.af.sum[k]
    surv.af[2*y,k] <- surv.af[2*y,k-1]*s.af.win[k-1]*s.af.sum[k]
    
    # model prob recovery = prob survives until k, recovered in k 
    pr.jf.shot[2*y-1,k] <- surv.jf[2*y-1,k]*f.af[k]
    pr.jf.shot[2*y,k]   <- surv.jf[2*y,k]*f.af[k]
    pr.af.shot[2*y-1,k] <- surv.jf[2*y-1,k]*f.af[k]
    pr.af.shot[2*y,k]   <- surv.jf[2*y,k]*f.af[k]
    } #k
    } # end second y loop
    
    # third y loop fills in remainder of dmulti matrix
    for (y in 1:yrs){  
    # Left of main diagonal (prob of encounter prior to banding = 0)
    for (l in 1:(y-1)){
    pr.jf.shot[2*y-1,l] <- 0
    pr.jf.shot[2*y,l] <- 0
    pr.af.shot[2*y-1,l] <- 0
    pr.af.shot[2*y,l] <- 0
    } #l
    } #y
    
    # Last column: probability of non-recovery
    for (t in 1:(2*yrs)){
    pr.jf.shot[t,(yrs+1)] <- 1-sum(pr.jf.shot[t,1:yrs])
    pr.af.shot[t,(yrs+1)] <- 1-sum(pr.af.shot[t,1:yrs])
    } #t
    }  # end bugs model  
    ",fill=TRUE)
sink()

# Bundle data 
bugs.data <- list(yrs=20, jf.obs.Aug=trunc(Lincoln.out[,1]/1000), af.obs.Aug=trunc(Lincoln.out[,5]/1000),
                  f.obs.Feb=trunc(Lincoln.out[,9]/1000), marr.juvF.shot=marr.juvF.shot, marr.ad_unkF.shot=marr.ad_unkF.shot)

# Initial values (starting chains in target range for each parameter)  
inits <- function() {list(s.af.sum.mu = runif(1,0,3), s.af.win.mu = runif(1,0,3), s.jf.win.mu = runif(1,0,3), 
                          f.af.mu = runif(1,-4,-1), f.jf.mu = runif(1,-4,-1), 
                          F.jf.mu = runif(1,0.5,2), N.f.Feb=c(rpois(1,958.3),rep(NA,19)))}

# Parameters monitored   
parameters <- c("N.f.Feb","N.af.Aug","N.jf.Aug",
                "s.af.sum.mu","s.af.win.mu","s.jf.win.mu",
                "s.af.sum.sd","s.af.win.sd","s.jf.win.sd",
                "f.af.mu","f.jf.mu","f.af.sd","f.jf.sd","F.jf.mu","F.jf.sd",
                "s.af.sum","s.af.win",,"s.jf.sum","s.jf.win",
                "f.af","f.jf","F.jf")

# MCMC settings (just trying to get this thing to go)
ni <- 100
nt <- 1
nb <- 0
nc <- 3

# Call WinBUGS from R (bugs model gives undefined real result, jags hangs on rel.jf[1])
ABDU.IPM2.out <- bugs(bugs.data, inits, parameters, "ABDU.IPM2.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,debug=TRUE,bugs.directory=bugs.dir,working.directory=getwd())
#print(ABDU.IPM2.out, digits=4)

ABDU.IPM2.out <- jags(bugs.data, inits, parameters, "ABDU.IPM2.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,working.directory=getwd())
#print(ABDU.IPM2.out, digits=4)


########################
# 3c training wheel version, S, f and associated variances from separate analysis of recovery data
# juv S during winter only
########################
Brownie <- read.csv("ABDU_surv_f.csv",header=TRUE)
head(Brownie)
tail(Lincoln.out)

sink("ABDU.IPM3.bug")
cat("
    model {
    
    #---------------------------------------------
    # 1. Priors and constraints
    #---------------------------------------------
    
    # observed initial Lincoln estimates/1000 (N) for combined juv + adult females (af) and males (am)
    N.f.Feb[1] ~ dpois(958.3)    
    N.m.Feb[1] ~ dpois(1033.9)
    
    F.jf.mu ~ dunif(0.5,2.5)      # prior for mean fecundity rates (F sex spec because fall sex ratio unbalanced)
    F.jm.mu ~ dunif(0.5,2.5)      
    F.jf.sd ~ dunif(0,2)                            
    F.jm.sd ~ dunif(0,2)    
    F.jf.tau <- pow(F.jf.sd,-2)    
    F.jm.tau <- pow(F.jm.sd,-2)    
    
    #---------------------------------------------
    # 2. Generate annual parameter estimates (S, f, F)
    #---------------------------------------------
    
    for (y in 1:yrs){
    s.af.sum[y] ~ dnorm(s.af.sum.mu[y],s.af.sum.tau[y])
    s.af.win[y] ~ dnorm(s.af.win.mu[y],s.af.win.tau[y])
    s.jf.win[y] ~ dnorm(s.jf.win.mu[y],s.jf.win.tau[y])
    s.am.sum[y] ~ dnorm(s.am.sum.mu[y],s.am.sum.tau[y])
    s.am.win[y] ~ dnorm(s.am.win.mu[y],s.am.win.tau[y])
    s.jm.win[y] ~ dnorm(s.jm.win.mu[y],s.jm.win.tau[y])
    Fecun.jf[y] ~ dnorm(F.jf.mu,F.jf.tau) 
    F.jf[y] <-max(0,Fecun.jf[y])                # avoiding I(0,) format
    Fecun.jm[y] ~ dnorm(F.jm.mu,F.jm.tau) 
    F.jm[y] <-max(0,Fecun.jm[y])
    }
    
    #-------------------------------------------
    # 3a. Likelihood for population components from Lincoln estimates
    #-------------------------------------------
    # millions of birds, so ignoring demographic stochasticity
    # estimate fall population sizes in year 1
    N.af.Aug[1] <- N.f.Feb[1] * s.af.sum[1]  
    N.am.Aug[1] <- N.m.Feb[1] * s.am.sum[1]  
    N.jf.Aug[1] <- N.af.Aug[1] * F.jf[1]
    N.jm.Aug[1] <- N.af.Aug[1] * F.jm[1]
    
    for (t in 2:yrs){
    N.f.Feb[t] <- (N.af.Aug[t-1] * s.af.win[t-1]) + (N.jf.Aug[t-1] * s.jf.win[t-1])
    N.m.Feb[t] <- (N.am.Aug[t-1] * s.am.win[t-1]) + (N.jm.Aug[t-1] * s.jm.win[t-1])
    N.af.Aug[t] <- N.f.Feb[t] * s.af.sum[t]  
    N.am.Aug[t] <- N.m.Feb[t] * s.am.sum[t]  
    N.jf.Aug[t] <- N.af.Aug[t] * F.jf[t]
    N.jm.Aug[t] <- N.af.Aug[t] * F.jm[t]
    }
    
    # Observation process, population counts
    for (t in 1:yrs) {
    jf.obs.Aug[t] ~ dpois(N.jf.Aug[t]) 
    jm.obs.Aug[t] ~ dpois(N.jm.Aug[t]) 
    af.obs.Aug[t] ~ dpois(N.af.Aug[t]) 
    am.obs.Aug[t] ~ dpois(N.am.Aug[t]) 
    f.obs.Feb[t] ~ dpois(N.f.Feb[t]) 
    m.obs.Feb[t] ~ dpois(N.m.Feb[t]) 
    }
    
    # Derived parameters, lambda
    for (t in 1:(yrs-1)) {
    lambda.f.Feb[t] <- N.f.Feb[t+1]/N.f.Feb[t] 
    lambda.m.Feb[t] <- N.m.Feb[t+1]/N.m.Feb[t] 
    lambda.af.Aug[t] <- N.af.Aug[t+1]/N.af.Aug[t] 
    lambda.am.Aug[t] <- N.am.Aug[t+1]/N.am.Aug[t] 
    }
    
    }  # end bugs model  
    ",fill=TRUE)
sink()

# Bundle data 
bugs.data <- list(yrs=20, jf.obs.Aug=trunc(Lincoln.out[,1]/1000), jm.obs.Aug=trunc(Lincoln.out[,3]/1000), af.obs.Aug=trunc(Lincoln.out[,5]/1000),
                  am.obs.Aug=trunc(Lincoln.out[,7]/1000), f.obs.Feb=trunc(Lincoln.out[,9]/1000), m.obs.Feb=trunc(Lincoln.out[,11]/1000),
                  s.jf.win.mu=Brownie$s.jf.win.mu, s.jf.win.tau=Brownie$s.jf.win.tau, s.jm.win.mu=Brownie$s.jm.win.mu, s.jm.win.tau=Brownie$s.jm.win.tau,
                  s.af.win.mu=Brownie$s.af.win.mu, s.af.win.tau=Brownie$s.af.win.tau, s.am.win.mu=Brownie$s.am.win.mu, s.am.win.tau=Brownie$s.am.win.tau,
                  s.af.sum.mu=Brownie$s.af.sum.mu, s.af.sum.tau=Brownie$s.af.sum.tau, s.am.sum.mu=Brownie$s.am.sum.mu, s.am.sum.tau=Brownie$s.am.sum.tau)

# Initial values (starting chains in target range) 
inits <- function() {list(F.jf.mu=runif(1,0.5,2.5), F.jm.mu = runif(1,0.5,2.5),      
                          F.jf.sd = runif(1,0,2), F.jm.sd = runif(1,0,2),    
                          N.f.Feb=c(rpois(1,958.3),rep(NA,19)),N.m.Feb=c(rpois(1,1033.9),rep(NA,19)))}

# Parameters monitored   
parameters <- c("N.f.Feb","N.m.Feb","N.af.Aug","N.am.Aug","N.jf.Aug","N.jm.Aug",
                "F.jf.mu","F.jm.mu","F.jf.sd","F.jm.sd",
                "s.af.sum","s.af.win","s.am.sum","s.am.win","s.jf.win","s.jm.win","F.jf","F.jm",
                "lambda.f.Feb","lambda.m.Feb","lambda.af.Aug","lambda.am.Aug")

# MCMC settings (just trying to get this thing to go)
ni <- 150000
nt <- 1
nb <- 50000
nc <- 3

# Call WinBUGS from R (bugs model gives undefined real result, jags hangs on rel.jf[1])
#ABDU.IPM3.out <- bugs(bugs.data, inits, parameters, "ABDU.IPM3.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,debug=TRUE,bugs.directory=bugs.dir,working.directory=getwd())
#print(ABDU.IPM.out, digits=4)

ABDU.IPM3.out <- jags(bugs.data, inits, parameters, "ABDU.IPM3.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,working.directory=getwd())
print(ABDU.IPM3.out, digits=4)
