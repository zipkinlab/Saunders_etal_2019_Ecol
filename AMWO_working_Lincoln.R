
# Make sure AMWO marray and AMWO rel are loaded into workspace
library(jagsUI)

# Load WING and HARVEST data #
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
clean[,9] <- rowSums(clean[,3:8])     # NOTE: Adjust if don't want all unknown groups

colnames(clean)<-c("Year","Pop","AdF","AdM","AdU","JuvF","JuvM","JuvU","Total")

wings<-matrix(NA,nrow=51,ncol=2)
wings[,1] <- clean$Total[clean$Pop=="E"]   # column 1 is eastern (matches p indexing below)
wings[,2] <-clean$Total[clean$Pop=="C"]    # column 2 is central

dim(wings)

# add 2 more rows for missing 2014-2015 years as placeholder for now
V1 <- c(rep(NA, 2))
V2 <- c(rep(NA, 2))
wings.add <- matrix(c(V1, V2), nrow=length(V1))
wings <- rbind(wings, wings.add)

# Now need to make wings.age and wings.sex datasets
clean[,10] <- rowSums(clean[,3:5])   # adult wings column; adjust categories included if needed
clean[,11] <- rowSums(clean[,6:8])   # juvenile wings column; adjust categories included if needed

colnames(clean)<-c("Year","Pop","AdF","AdM","AdU","JuvF","JuvM","JuvU","Total", "Adults", "Juvs")

wings.age<-array(NA, dim=c(51,2,2))              # dims: t, ages (1=juv; 2=adult), p (1=eastern; 2=central)
wings.age[,1,1]<-clean$Juvs[clean$Pop=="E"]      #juvenile eastern
wings.age[,2,1]<-clean$Adults[clean$Pop=="E"]    #adult eastern
wings.age[,1,2]<-clean$Juvs[clean$Pop=="C"]      #juvenile central
wings.age[,2,2]<-clean$Adults[clean$Pop=="C"]    # adult central

dim(wings.age)    #51 by 2 by 2 (need to add more years)

# add 2 more rows for missing 2014-2015 years as placeholder for now
wings.new <- array(NA, dim=c(53,2,2))
wings.new[-52:-53,,] <- wings.age
wings.age <- wings.new

dim(wings.age)

# wings.sex dataset (adults only)
wings.sex<-array(NA, dim=c(51,2,2))              # dims: t, sex (1=male; 2=female), p (1=eastern; 2=central)
wings.sex[,1,1]<-clean$AdM[clean$Pop=="E"]      # male eastern
wings.sex[,2,1]<-clean$AdF[clean$Pop=="E"]    # female eastern
wings.sex[,1,2]<-clean$AdM[clean$Pop=="C"]      # male central
wings.sex[,2,2]<-clean$AdF[clean$Pop=="C"]    # female central

dim(wings.sex)  #51 by 2 by 2 (need to add more years)

# add 2 more rows for missing years as placeholder for now
wings.new <- array(NA, dim=c(53,2,2))
wings.new[-52:-53,,] <- wings.sex
wings.sex <- wings.new

dim(wings.sex)

# Bring in total harvest data
#add total harvest column to clean. 
############
##NOTE: this is Harvest_Est column using harvest estimates from Todd's model!#####
############
clean[,12] <- harvest$Harvest_Est
colnames(clean)<-c("Year","Pop","AdF","AdM","AdU","JuvF","JuvM","JuvU","Total","Adults","Juvs","Harvest")

H.total<-matrix(NA, nrow=51, ncol=2)          #dims: t by p (1=eastern; 2=central)
H.total[,1]<-clean$Harvest[clean$Pop=="E"]    # harvest in eastern pop
H.total[,2]<-clean$Harvest[clean$Pop=="C"]    # harvest in central pop

dim(H.total)      #51 by 2

# add 2 more rows for missing 2014-2015 years as placeholder for now
V1 <- c(rep(NA, 2))
V2 <- c(rep(NA, 2))
H.add <- matrix(c(V1, V2), nrow=length(V1))
H.total <- rbind(H.total, H.add)

dim(H.total)  #53 by 2

## Ones trick vector ##################################
ones <- array(1, dim = c(53, 3, 2)) 


#------------#
#-BUGS Model-#
#------------#

sink("Lincoln.Simple2.6April.jags")
cat("
    model {
    
    # Priors and constraints for population means and variances
    
    # prior for initial pop sizes in spring               # informed prior based on Lincoln estimates
    n[1,1,3,1] ~ dunif(1000000, 3000000)                  #dnorm(700000,1E-10)I(0,)     # initial size for female eastern
    n[1,1,3,2] ~ dunif(1000000, 3000000)                  #dnorm(1600000,5E-12)I(0,)    # initial size for female central
    n[1,1,2,1] ~ dunif(1000000, 3000000)                  #dnorm(530000,1E-10)I(0,)     # initial size for male eastern
    n[1,1,2,2] ~ dunif(1000000, 3000000)                  #dnorm(615000,1E-10)I(0,)     # initial size for male central
    n[1,1,1,1] <- n[1,1,3,1] * F[1,1]                        # initial size for juv eastern
    n[1,1,1,2] <- n[1,1,3,2] * F[1,2]                        # initial size for juv central
    
    # round initial sizes
    N[1,1,3,1] <- round(n[1,1,3,1])
    N[1,1,3,2] <- round(n[1,1,3,2])
    N[1,1,2,1] <- round(n[1,1,2,1])
    N[1,1,2,2] <- round(n[1,1,2,2])
    N[1,1,1,1] <- round(n[1,1,1,1])
    N[1,1,1,2] <- round(n[1,1,1,2])
    
    pi.sex[1] <- pi[2]/(pi[2]+pi[3])           # male to female proportion
    pi.sex[2] <- 1-pi.sex[1]                   # female to male proportion
    
    #for (t in 1:yrs){
    #report[t] ~ dunif(0.15, 0.85)                       # prior for reporting rate (for now)
    #} #t
    
    for (p in 1:2){                            # 1 Eastern, 2 Central
    
    # priors for fecundities
    F.x[p] ~ dunif(0, 4)                       # logical bounds on fecundity (F) with 4 egg clutch 
    F.sd[p] ~ dunif(0.01, 1) 
    F.tau[p] <- pow(F.sd[p], -2)
    
    for (t in 1:yrs){
    F[t,p] ~ dnorm(F.x[p], F.tau[p])T(0,)                        # Fecundity cannot be <0 or >4
    } #t
    
    for (c in 1:3){                       
    sa.x[c,p] ~ dunif(0,1) 
    sa.mu[c,p] <- logit(sa.x[c,p])        
    ss.x[c,p] ~ dunif(0,1) 
    ss.mu[c,p] <- logit(ss.x[c,p])
    f.x[c,p] ~ dunif(0,0.1) 
    f.mu[c,p] <- logit(f.x[c,p])               
    
    sa.sd[c,p] ~ dunif(0.05,2)                 # Priors for SDs of survival and recovery rates
    ss.sd[c,p] ~ dunif(0.05,2) 
    f.sd[c,p] ~ dunif(0.05,2) 
    sa.tau[c,p] <- pow(sa.sd[c,p],-2)          # Express as precision
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
    
    } # close t
    
    ## Balance equations, assuming estimating N in both seasons
    # Initial pop sizes in summer [t, s, c, p]
    N[1,2,c,p] <- trunc(N[1,1,c,p]*ss[1,c,p])         # pop size in late summer is function of summer survival 
    } #c
    
    for (t in 2:yrs){
    N[t,1,1,p] <- trunc(N[t,1,3,p] * F[t,p])          # pop size of juv in spring = adF * mean fecundity
    N[t,2,1,p] <- trunc(N[t,1,1,p]*ss[t,1,p])         # pop size of juv in late summer
    sw[t,1,p] <- sa[t-1,1,p]/ss[t,1,p]         # derived winter survival for juveniles
    
    for (c in 2:3){
    sw[t,c,p] <- sa[t-1,c,p]/ss[t,c,p]                                            # derived winter survival for adults
    N[t,1,c,p] <- trunc(N[t-1,2,c,p]*sw[t,c,p] + pi.sex[c-1]*N[t-1,2,1,p]*sw[t,1,p])     # pop size of adults in spring
    N[t,2,c,p] <- trunc(N[t,1,c,p]*ss[t,c,p])                                            # pop size of adults in late summer 
    } #c
    } #t 
    
    for (t in 1:yrs){
    # Harvest estimates by age-sex class are function of harvest rate and total pop size
    for (c in 1:3){

    H[t,c,p] ~ dbin(h[t,c,p], N[t,2,c,p])
    
    h[t,c,p] <- f[t,c,p]/0.506                   # harvest rate (h) is recovery rate divided by reporting rate (p)                                                      
    } #c                                              # note that report is likely going to be multipled by vector of proportions of 1800 bands used each year
    } #t           
    } #p
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data

pi[1] <- 0.46
pi[2] <- 0.21
pi[3] <- 0.33

H <- array(NA, dim = c(53,3,2))
for(t in 1:53){
  for(c in 1:3){
    for(p in 1:2){
      H[t,c,p] <- round(pi[c] * H.total[t,p])
    }
  }
}

H[52:53,1,] <- as.integer(round(mean(H[,1,], na.rm = TRUE)))
H[52:53,2,] <- as.integer(round(mean(H[,2,], na.rm = TRUE)))
H[52:53,3,] <- as.integer(round(mean(H[,3,], na.rm = TRUE)))



bugs.data <- list(yrs=53, H=H, pi=pi)  

# inits

#N.inits <- array(NA, dim=c(53, 2, 3, 2))
#for (t in 2:53){
#for (c in 1:3){
#for (p in 1:2){
#for (s in 1:2){
#N.inits[t,s,c,p] <- round(runif(1, 100000, 1000000))
#}}}}

H.inits <- array(NA, dim=c(53, 3, 2))   #need this to be by t as well?
for (t in 1:53){
  H.inits[t,1,1] <- round(rnorm(1,96000, 10000))    #juv eastern
  H.inits[t,1,2] <- round(rnorm(1,177000, 20000))   #juv central
  H.inits[t,2,1] <- round(rnorm(1,40000, 5000))    # male eastern
  H.inits[t,2,2] <- round(rnorm(1,73000, 10000))   # male central
  H.inits[t,3,1] <- round(rnorm(1,60000, 10000))   # female eastern
  H.inits[t,3,2] <- round(rnorm(1,116000, 20000))   #female central
}

f.x.inits <- matrix(NA, nrow=3, ncol=2)
for (c in 1:3){
  for (p in 1:2){
    f.x.inits[c,p] <- runif(1, 0, 0.05)
  }}

f.sd.inits <- matrix(NA, nrow=3, ncol=2)
for (c in 1:3){
  for (p in 1:2){
    f.sd.inits[c,p] <- runif(1, 0.05, 1.5)
  }}

sa.x.inits <- matrix(NA, nrow=3, ncol=2)                                                 
fill <- runif(1,0.2,0.9)                                     
for (p in 1:2){
  for (c in 1:3){
    sa.x.inits[c,p] <- fill
  }}

ss.x.inits <- matrix(NA, nrow=3, ncol=2)
fill <- runif(1,0.1,0.6)
for (p in 1:2){
  for (c in 1:3){
    ss.x.inits[c,p] <- fill
  }}

report.inits <- rep(NA, 53)
for (t in 1:53){
  report.inits[t] <- runif(1, 0.4, 0.85)
}

n.inits <- array(NA, dim=c(1,1,3,2))
n.inits[1,1,3,1] <- 2000000
n.inits[1,1,3,2] <- 1600000
n.inits[1,1,2,1] <- 1200000
n.inits[1,1,2,2] <- 1700000

F.x.inits <- rep(NA, 2)
for (p in 1:2){
  F.x.inits[p] <- 2.5
}

inits <- function(){list(#H=H.inits,
  #pi.sex=pi.sex.inits, 
  #sa.x=sa.x.inits, 
  #ss.x=ss.x.inits,  
  #f.sd=f.sd.inits,
  #f.x=f.x.inits, 
  #F.x=F.x.inits,
  n=n.inits
)} #,report=report.inits)}                 

# Parameters monitored
parameters <- c("sa.x", "ss.x", "f.x", "sa.sd", "ss.sd", "f.sd", "sa", "f", "N")  # add to this before running!

# MCMC settings
ni <- 1000
nt <- 1
nb <- 100
nc <- 1

# Run JAGS
AMWO.lincoln.jags <- jagsUI(bugs.data, inits=inits, parameters, "Lincoln.Simple2.6April.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE, store.data=TRUE) #changed from inits=NULL
