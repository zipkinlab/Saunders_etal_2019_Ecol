########################################################################
# Combined Lincoln estimation-Brownie recovery model for AmWO
########################################################################

#-----------------------------------------------------------------------------------------------
# Indexing: 
#t for cohort year
#k for recovery year (if NE to cohort year)
#s for season [1=Apr-Jun, 2=Jul-Sep]
#c for cohort [1=juvenile, 2=adult male, 3=adult female]
#p for population [1=eastern, 2=central]
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# Things we still have to do:
# 1. Get wing and harvest data read in and manipulated into correct format to use 3/23/17: done
# 2. Figure out if we can use pi.age to inform fecundity? 3/23/17: worry about this later
# 3. Figure out initial pop sizes to use  3/23/17: used Todd's emailed pop sizes
# 4. Figure out what we want reporting rate (p) to be a product of: vector of length yrs inc.
# proportion of 1800 bands and band type effect?  3/23/17: using semi-informative prior for now
# 5. Need missing harvest/wing data for 1997, 1998, 2014 and 2015  3/23/17: putting in NAs for now
#-----------------------------------------------------------------------------------------------

# Make sure AMWO marray and AMWO rel are loaded into workspace

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

#add total harvest column to clean
clean[,12] <- harvest$Harvest
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

#------------#
#-BUGS Model-#
#------------#

sink("Lincoln.Brownie.23March.jags")
cat("
    model {
  
    # Priors and constraints for population means and variances

    # prior for initial pop sizes in spring             # informed prior based on Lincoln estimates
    n[1,1,3,1] ~ dnorm(600000,1E-10)                    # initial size for female eastern
    n[1,1,3,2] ~ dnorm(720000,1E-10)                    # initial size for female central
    n[1,1,2,1] ~ dnorm(530000,1E-10)                    # initial size for male eastern
    n[1,1,2,2] ~ dnorm(615000,1E-10)                    # initial size for male central
    n[1,1,1,1] <- n[1,1,3,1] * F[1,1]                   # initial size for juv eastern
    n[1,1,1,2] <- n[1,1,3,2] * F[1,2]                   # initial size for juv central

    # round initial sizes
    N[1,1,3,1] <- round(n[1,1,3,1])
    N[1,1,3,2] <- round(n[1,1,3,2])
    N[1,1,2,1] <- round(n[1,1,2,1])
    N[1,1,2,2] <- round(n[1,1,2,2])
    N[1,1,1,1] <- round(n[1,1,1,1])
    N[1,1,1,2] <- round(n[1,1,1,2])

    for (t in 1:yrs){
    report[t] ~ dunif(0.15, 0.85)                       # prior for reporting rate (for now)
    } #t

    for (p in 1:2){                            # 1 Eastern, 2 Central

    # priors for fecundities
    F.x[p] ~ dunif(0, 4)                       # logical bounds on fecundity (F) with 4 egg clutch 
    F.sd[p] ~ dunif(0.01, 1) 
    F.tau[p] <- pow(F.sd[p], -2)

    for (t in 1:yrs){
    F[t,p] ~ dnorm(F.x[p], F.tau[p])                       # Fecundity cannot be <0 or >4
    } #t

    for (c in 1:3){                       
    sa.x[c,p] ~ dunif(0,1) 
    sa.mu[c,p] <- logit(sa.x[c,p])        
    ss.x[c,p] ~ dunif(0,1) 
    ss.mu[c,p] <- logit(ss.x[c,p])
    f.x[c,p] ~ dunif(0,1) 
    f.mu[c,p] <- logit(f.x[c,p])               # note: may need a prior for pi.sex (see below for where pi.sex created)

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
    N[1,2,c,p] <- N[1,1,c,p]*ss[1,c,p]         # pop size in late summer is function of summer survival 
    } #c

    for (t in 2:yrs){
    N[t,1,1,p] <- N[t,1,3,p] * F[t,p]          # pop size of juv in spring = adF * mean fecundity
    N[t,2,1,p] <- N[t,1,1,p]*ss[t,1,p]         # pop size of juv in late summer
    sw[t,1,p] <- sa[t-1,1,p]/ss[t,1,p]         # derived winter survival for juveniles

    for (c in 2:3){
    sw[t,c,p] <- sa[t-1,c,p]/ss[t,c,p]                                                # derived winter survival for adults
    N[t,1,c,p] <- N[t-1,2,c,p]*sw[t,c,p] + pi.sex[t,c-1,p]*N[t-1,2,1,p]*sw[t,1,p]     # pop size of adults in spring
    N[t,2,c,p] <- N[t,1,c,p]*ss[t,c,p]                                                # pop size of adults in late summer 
    } #c
    } #t 

    # Observation process, recoveries and harvest data
    # Recoveries
    # Note: releases MUST be provided as data; we have saved as relAMWO
    for (c in 1:3){
    for (t in 1:yrs){
    marrayAMWO[t,1:(yrs+1),1,c,p] ~ dmulti(pr[t,,1,c,p], rel[t,1,c,p])        # Apr-Jun releases     
    marrayAMWO[t,1:(yrs+1),2,c,p] ~ dmulti(pr[t,,2,c,p], rel[t,2,c,p])        # Jul-Sep releases
    
    ## Define cell probabilities of m-arrays
    # Main diagonal--direct recovery in season [t]
    pr[t,t,1,c,p] <- ss[t,c,p] * f[t,c,p]                  # spring banded birds must survive to start of hunting season
    pr[t,t,2,c,p] <- f[t,c,p]                              # survival assumed 1 if banded Jul-Sep
    
    # monitor cumulative survival to start of next hunting season (previously surv; used in subsequent diagonals)
    cumS[t,t,1,c,p] <- ss[t,c,p] * sa[t,c,p]
    cumS[t,t,2,c,p] <- sa[t,c,p]
    } # t
    } # c
    
    for (t in 1:yrs){
    # Above main diagonal--indirect recovery in season [k > t]
    # All birds are adults now, no differences between Apr-Jun and Jul-Sep either
    for (k in (t+1):yrs){                                                                       # k loop to represent next year (above main diagonal)
    # recoveries
    pr[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi.sex[t,1,p]*f[k,2,p] + pi.sex[t,2,p]*f[k,3,p])     # juvs as mixture of AdM & AdF
    pr[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi.sex[t,1,p]*f[k,2,p] + pi.sex[t,2,p]*f[k,3,p])    # pi.sex[,1,] is proportion male; pi.sex[,2,] is proportion female. See below       
    pr[t,k,1,2,p] <- cumS[t,k-1,1,2,p] * f[k,2,p] 
    pr[t,k,2,2,p] <- cumS[t,k-1,2,2,p] * f[k,2,p] 
    pr[t,k,1,3,p] <- cumS[t,k-1,1,3,p] * f[k,3,p] 
    pr[t,k,2,3,p] <- cumS[t,k-1,2,3,p] * f[k,3,p] 
    # monitor cumulative survival to start of next hunting period
    cumS[t,k,1,1,p] <- cumS[t,k-1,1,1,p] * (pi.sex[t,1,p]*sa[k,2,p] + pi.sex[t,2,p]*sa[k,3,p]) # juvs as mixture AdM & AdF
    cumS[t,k,2,1,p] <- cumS[t,k-1,2,1,p] * (pi.sex[t,1,p]*sa[k,2,p] + pi.sex[t,2,p]*sa[k,3,p]) 
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

    # Summarize known wing samples by age and sex
    for (t in 1:yrs){
    wings.age[t,1,p] ~ dbin(pi.age[t,1,p], wings[t,p])          # pi.age is the proportion of juveniles to adults    
    wings.sex[t,,p] ~ dmulti(pi.sex[t,,p], wings.age[t,2,p])    # this needs to stay dmulti because of pi.sex on left-side of equations below; can't put (1-pi.sex[t,1,p]) as a derived quantity    

    # Bring in total harvest data
    pi.age[t,1,p] <- H[t,1,p]/H.total[t,p]                          # pi is on the left-hand side (instead of H) because we need to use H on the left-hand side below
    pi.sex[t,1,p] <- H[t,2,p]/((1-pi.age[t,1,p]) * H.total[t,p])    # can we use pi.age to inform fecundity??
    pi.sex[t,2,p] <- H[t,3,p]/((1-pi.age[t,1,p]) * H.total[t,p]) 

    # Harvest estimates by age-sex class are function of harvest rate and total pop size
    for (c in 1:3){
    H[t,c,p] ~ dbin(h[t,c,p], N[t,2,c,p]) 
    #H[t,c,p] <- h[t,c,p]*N[t,2,c,p]                  # use this instead of dbin above?
    h[t,c,p] <- f[t,c,p]/report[t]                    # harvest rate (h) is recovery rate divided by reporting rate (p)                                                      
    } #c                                              # note that p is likely going to be multipled by vector of proportions of 1800 bands used each year
    } #t           
    } #p
    } # end bugs model
    ",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(yrs=dim(marrayAMWO)[1], marrayAMWO=marrayAMWO, rel=relAMWO, wings.age=wings.age, wings.sex=wings.sex, H.total=H.total, wings=wings)  

### inits not needed if they are mirror image of priors
### but "moderately informed inits" can become essential to get models running, so good to have    

N.inits <- array(NA, dim=c(53, 2, 3, 2))
for (t in 1:53){
  for (c in 1:3){
    for (p in 1:2){
      for (s in 1:2){
      N.inits[t,s,c,p] <- round(runif(1, 100000, 1000000))
    }}}}

H.inits <- array(NA, dim=c(53, 3, 2))
for (t in 1:53){
  for (c in 1:3){
    for (p in 1:2){
      H.inits[t,c,p] <- round(runif(1,50000, 600000))
    }}}

inits <- function(){list(N=N.inits, H=H.inits)}

#sa.x.inits <- matrix(NA, nrow=3, ncol=2)                                                 
#fill <- runif(1,0,1)                                     
#for (p in 1:2){
  #for (c in 1:3){
    #sa.x.inits[c,p] <- fill
  #}}

#ss.x.inits <- matrix(NA, nrow=3, ncol=2)
#fill <- runif(1,0,1)
#for (p in 1:2){
  #for (c in 1:3){
    #ss.x.inits[c,p] <- fill
  #}}

#f.x.inits <- matrix(NA, nrow=3, ncol=2)
#fill <- runif(1,0,1)
#for (p in 1:2){
  #for (c in 1:3){
    #f.x.inits[c,p] <- fill
  #}}

#sa.sd.inits <- matrix(NA, nrow=3, ncol=2)
#fill <- runif(1,0.05,2)
#for (p in 1:2){
  #for (c in 1:3){
    #sa.sd.inits[c,p] <- fill
  #}}

#ss.sd.inits <- matrix(NA, nrow=3, ncol=2)
#fill <- runif(1,0.05,2)
#for (p in 1:2){
  #for (c in 1:3){
    #ss.sd.inits[c,p] <- fill
  #}}

#f.sd.inits <- matrix(NA, nrow=3, ncol=2)
#fill <- runif(1,0.05,2)
#for (p in 1:2){
  #for (c in 1:3){
    #f.sd.inits[c,p] <- fill
  #}}

# Parameters monitored
parameters <- c("sa.x", "ss.x", "f.x", "sa.sd", "ss.sd", "f.sd", "sa", "f", "pi.sex", "pi.age")  # add to this before running!

# MCMC settings
ni <- 1000
nt <- 1
nb <- 100
nc <- 1

# Run JAGS
AMWO.combo.jags <- jagsUI(bugs.data, inits=inits, parameters, "Lincoln.Brownie.23March.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE) #changed from inits=NULL


