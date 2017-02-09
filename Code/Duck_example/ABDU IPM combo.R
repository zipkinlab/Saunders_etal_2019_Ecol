#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("C:/Users/farrm/Documents/GitHub/timberdoodle/Code/Duck_example")

#----------------#
#-Load libraries-#
#----------------#

library(jagsUI)

#-----------#
#-Load data-#
#-----------#

harvest <- read.csv("ABDU harvest and wing data.csv",header=TRUE)

banding <- read.csv("ABDU shot only.csv",header=FALSE)

#-----------------------------------------#
#-Create matrices for Lincoln estimation-#
#-----------------------------------------#

#Matrix to hold US harvest data from 1969 to 1988
US.harvest <- matrix(0,nrow=20,ncol=13,byrow=TRUE)

#Matrix to hold Can harvest data from 1969 to 1988
Can.harvest <- matrix(0,nrow=20,ncol=13,byrow=TRUE)

#Matrix to hold US wing data from 1969 to 1988
US.wings <- matrix(0,nrow=20,ncol=8,byrow=TRUE)

#Matrix to hold Can wing data from 1969 to 1988
Can.wings <- matrix(0,nrow=20,ncol=8,byrow=TRUE)

for (i in 1:20){
  for (k in 1:13){
    
    #Rows 1:20 US; Col 1:13 harvest data
    US.harvest[i,k] <- harvest[i,k]
    
    #Rows 21:40 Can; Col 1:13 harvest data
    Can.harvest[i,k] <- harvest[i+20,k]
  }
  for (j in 1:6){
    #Rows 1:20 US; Col 14:20 wing data
    US.wings[i,j] <- harvest[i,j+13]
    
    #Rows 21:40 Can; Col 14:20 wing data
    Can.wings[i,j] <- harvest[i+20,j+13]
  } #End j
  
  #Summarize known adults (including unknown sex)
  US.wings[i,7] <- sum(US.wings[i,1:3])
  
  #Summarize known juvs (including unknown sex)
  US.wings[i,8] <- sum(US.wings[i,4:6]) 
  
  #Summarize known adults (including unknown sex)
  Can.wings[i,7] <- sum(Can.wings[i,1:3])
  
  #Summarize known juvs (including unknown sex)
  Can.wings[i,8] <- sum(Can.wings[i,4:6]) 
} #End i

#Matrix to hold total bandings; populated in m-arrays
banded <- matrix(0,nrow=40,ncol=8,byrow=TRUE)

#-------------------------#
#-Create banding m-arrays-#
#-------------------------#

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

#------------#
#-BUGS Model-#
#------------#

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
    #SS: we will be estimating a combined 'correction factor' instead of values below
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
    
    #Calculate age and sex probability
    for (y in 1:yrs){
    US.wings[y,8] ~ dbin(pi.US.juv[y],US.known.age[y])
    US.wings[y,4] ~ dbin(pi.US.juvF[y],US.sexed.juv[y])
    US.wings[y,1] ~ dbin(pi.US.adF[y],US.sexed.ad[y])
    Can.wings[y,8] ~ dbin(pi.Can.juv[y],Can.known.age[y])
    Can.wings[y,4] ~ dbin(pi.Can.juvF[y],Can.sexed.juv[y])
    Can.wings[y,1] ~ dbin(pi.Can.adF[y],Can.sexed.ad[y])
    }
    
    #Generate component estimates of harvest from age and sex probability
    for (y in 1:yrs){
    Harv.JuvF.Can[y] <- pi.Can.juv[y]*pi.Can.juvF[y]*Can.Harvest[y,12]
    Harv.JuvM.Can[y] <- pi.Can.juv[y]*(1-pi.Can.juvF[y])*Can.Harvest[y,12]
    Harv.AdF.Can[y] <- (1-pi.Can.juv[y])*pi.Can.adF[y]*Can.Harvest[y,12]
    Harv.AdM.Can[y] <- (1-pi.Can.juv[y])*(1-pi.Can.adF[y])*Can.Harvest[y,12]
    Harv.JuvF.US[y] <- pi.US.juv[y]*pi.US.juvF[y]*US.Harvest[y,12]*harv.adj
    Harv.JuvM.US[y] <- pi.US.juv[y]*(1-pi.US.juvF[y])*US.Harvest[y,12]*harv.adj
    Harv.AdF.US[y] <- (1-pi.US.juv[y])*pi.US.adF[y]*US.Harvest[y,12]*harv.adj
    Harv.AdM.US[y] <- (1-pi.US.juv[y])*(1-pi.US.adF[y])*US.Harvest[y,12]*harv.adj
    

    H.juvF.fall[y] <- round(Harv.JuvF.Can[y] + Harv.JuvF.US[y])
    H.adF.fall[y] <- round(Harv.AdF.Can[y] + Harv.AdF.US[y])
    H.juvM.fall[y] <- round(Harv.JuvM.Can[y] + Harv.JuvM.US[y])
    H.adM.fall[y] <- round(Harv.AdM.Can[y] + Harv.AdM.US[y])
    H.combF.spring[y] <- round(Harv.AdF.Can[y] + Harv.AdF.US[y])
    H.combM.spring[y] <- round(Harv.AdM.Can[y] + Harv.AdM.US[y])
    }
    
    # Generate direct recovery rates
    for (y in 1:yrs){
    marr.juvF.shot.R[y*2,y] ~ dbin(f.pre.juvF[y],banded[y*2,1])
    marr.adF.shot.R[y*2,y] ~ dbin(f.pre.adF[y],banded[y*2,2])
    marr.combF.shot.R[y*2-1,y] ~ dbin(f.post.combF[y],banded[y*2-1,4])
    marr.juvM.shot.R[y*2,y] ~ dbin(f.pre.juvM[y],banded[y*2,5])
    marr.adM.shot.R[y*2,y] ~ dbin(f.pre.adM[y],banded[y*2,6])
    marr.combM.shot.R[y*2-1,y] ~ dbin(f.post.combM[y],banded[y*2-1,8])
    }
    
    # Generate Lincoln estimates and population structure as derived parameters 
    for (y in 1:yrs){
    
    h.pre.juvF[y] <- f.pre.juvF[y]/report
    h.pre.adF[y] <- f.pre.adF[y]/report
    h.post.combF[y] <- f.post.combF[y]/report
    h.pre.juvM[y] <- f.pre.juvM[y]/report
    h.pre.adM[y] <- f.pre.adM[y]/report
    h.post.combM[y] <- f.post.combM[y]/report

    H.juvF.fall[y] ~ dbin(h.pre.juvF[y], N.juvF.fall[y])
    H.adF.fall[y] ~ dbin(h.pre.adF[y], N.adF.fall[y])
    H.combF.spring[y] ~ dbin(h.post.combF[y], N.combF.spring[y])
    H.juvM.fall[y] ~ dbin(h.pre.juvM[y], N.juvM.fall[y])
    H.adM.fall[y] ~ dbin(h.pre.adM[y], N.adM.fall[y])
    H.combM.spring[y] ~ dbin(h.post.combM[y], N.combM.spring[y])

    #N.total.fall[y] <- N.juvF.fall[y] + N.juvM.fall[y] + N.adF.fall[y] + N.adM.fall[y]
    #N.total.spring[y] <- N.combF.spring[y] + N.combM.spring[y]
    #fecundity[y] <- 0.5*(N.juvF.fall[y]+N.juvM.fall[y])/N.adF.fall[y]
    #sex.ratio.juv[y] <- N.juvM.fall[y]/(N.juvF.fall[y]+N.juvM.fall[y])
    #sex.ratio.ad.fall[y] <- N.adM.fall[y]/(N.adF.fall[y]+N.adM.fall[y])
    #sex.ratio.comb.spring[y] <- N.combM.spring[y]/(N.combF.spring[y]+N.combM.spring[y])
    
    }

############################################################################################    
#---------------------------------------------
# 1. Priors and constraints
#---------------------------------------------

# observed initial Lincoln estimates/1000 (N) for combined juv + adult females (af) and males (am)
N.combF.spring[1] ~ dpois(958.3)    
N.combM.spring[1] ~ dpois(1033.9)

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

#DELETE?
#f.af.mu ~ dunif(0,0.3)                            
#lf.af.mu <- log(f.af.mu/(1-f.af.mu))    
#f.am.mu ~ dunif(0,0.3)                            
#lf.am.mu <- log(f.am.mu/(1-f.am.mu))    
#f.jf.mu ~ dunif(0,0.3)                            
#lf.jf.mu <- log(f.jf.mu/(1-f.jf.mu))    
#f.jm.mu ~ dunif(0,0.3)                            
#lf.jm.mu <- log(f.jm.mu/(1-f.jm.mu))    

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

#DELETE?
#f.af.sd ~ dunif(0.05,2)                            
#f.am.sd ~ dunif(0.05,2)                            
#f.jf.sd ~ dunif(0.05,2)                            
#f.jm.sd ~ dunif(0.05,2)                            

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

#DELETE?
#f.af.tau <- pow(f.af.sd,-2)                            
#f.am.tau <- pow(f.am.sd,-2)                            
#f.jf.tau <- pow(f.jf.sd,-2)                            
#f.jm.tau <- pow(f.jm.sd,-2)                            

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
  
#DELETE?
  #logit.f.af[y] ~ dnorm(lf.af.mu,f.af.tau)
  #logit(f.af[y]) <- logit.f.af[y] 
  #logit.f.am[y] ~ dnorm(lf.am.mu,f.am.tau)
  #logit(f.am[y]) <- logit.f.am[y] 
  #logit.f.jf[y] ~ dnorm(lf.jf.mu,f.jf.tau)
  #logit(f.jf[y]) <- logit.f.jf[y] 
  #logit.f.jm[y] ~ dnorm(lf.jm.mu,f.jm.tau)
  #logit(f.jm[y]) <- logit.f.jm[y] 
  

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
N.adF.spring[1] <- round(N.combF.spring[1] * (1-0.57))  #SS: what's this 0.57 and 0.42 exactly?
N.adM.spring[1] <- round(N.combM.spring[1] * (1-0.42))  
N.juvF.spring[1] <- round(N.combF.spring[1] * 0.57)
N.juvM.spring[1] <- round(N.combM.spring[1] * 0.42)


# estimate fall population sizes in year 1
N.adF.fall[1] <- round(N.adF.spring[1] * s.af.sum[1])  
N.adM.fall[1] <- round(N.adM.spring[1] * s.am.sum[1])  
N.juvF.fall[1] <- round(N.adF.fall[1] * F.jf[1])
N.juvM.fall[1] <- round(N.adF.fall[1] * F.jm[1])

for (t in 2:yrs){
N.adF.spring[t] <- N.adF.fall[t-1] * s.af.win[t-1]
N.adM.spring[t] <- N.adM.fall[t-1] * s.am.win[t-1]
N.juvF.spring[t] <- N.juvF.fall[t-1] * s.jf.win[t-1]
N.juvM.spring[t] <- N.juvM.fall[t-1] * s.jm.win[t-1]
# combined ad + juv population is what Lincoln (and BPOP) estimates
N.combF.spring[t] <- N.adF.spring[t] + N.juvF.spring[t]        
N.combM.spring[t] <- N.adM.spring[t] + N.juvM.spring[t]        
# juveniles (yearlings) graduate to adults in this step
N.adF.fall[t] <- N.adF.spring[t] * s.af.sum[t]  + N.juvF.spring[t] * s.jf.sum[t]
N.adM.fall[t] <- N.adM.spring[t] * s.am.sum[t]  + N.juvM.spring[t] * s.jm.sum[t]
N.juvF.fall[t] <- N.adF.fall[t] * F.jf[t]
N.juvM.fall[t] <- N.adF.fall[t] * F.jm[t]
}

# Derived parms, more efficient to estimate from saved posteriors unless I want intrinsic correlations
#  lambda.BPOP[t-1] <- (N.f.Feb[t]+N.m.Feb[t])/(N.f.Feb[t-1]+N.m.Feb[t-1])
#  lambda.Fall[t-1] <- (N.af.Aug[t]+N.am.Aug[t]+N.jf.Aug[t]+N.jm.Aug[t])/(N.af.Aug[t-1]+N.am.Aug[t-1]+N.jf.Aug[t-1]+N.jm.Aug[t-1])
#  sexratio.BPOP[t] <- (N.m.Feb[t])/(N.f.Feb[t]+N.m.Feb[t])

#-------------------------------------------
# 3a. Likelihood for band recovery data (S and f)
#-------------------------------------------

#SS: this is the Brownie recovery model as above, just incorporated into IPM

# Generate multinomial likelihoods for m-arrays (modified from Kery & Schaub)
# 2 rows for each year, Feb bandings, then Aug bandings
# also have m-arrays for other dead encounters and live recaps, but they add little precision, shot only here

for (t in 1:(2*yrs)){
  marr.juvF.shot[t,1:(yrs+1)] ~ dmulti(pr.jf.shot[t,], rel.jf[t])   #check/compare m-arrays
  marr.ad_unkF.shot[t,1:(yrs+1)] ~ dmulti(pr.af.shot[t,], rel.af[t])
  marr.juvM.shot[t,1:(yrs+1)] ~ dmulti(pr.jm.shot[t,], rel.jm[t])
  marr.ad_unkM.shot[t,1:(yrs+1)] ~ dmulti(pr.am.shot[t,], rel.am[t])
}

# Define cell probabilities of the m-arrays
for (y in 1:yrs){  

  #Direct recovery; main diagnol

#Cell prob for Juv Females Post-Harvest
pr.jf.shot[2*y-1,y] <- s.jf.sum[y] * f.pre.adF[y]

#Cell prob for Juv Females Pre-Harvest
pr.jf.shot[2*y,y] <- 1 * f.pre.juvF[y]              

#Cell prob for Adult Females Post-Harvest
pr.af.shot[2*y-1,y] <- s.af.sum[y] * f.pre.adF[y]

#Cell prob for Adult Females Pre-Harvest
pr.af.shot[2*y,y] <- 1 * f.pre.adF[y]

#Cell prob for Juv Males Post-Harvest
pr.jm.shot[2*y-1,y] <- s.jm.sum[y] * f.pre.adM[y]

#Cell prob for Juv Males Pre-Harvest
pr.jm.shot[2*y,y] <- 1 * f.pre.juvM[y]

#Cell prob for Adult Males Post-Harvest
pr.am.shot[2*y-1,y] <- s.am.sum[y] * f.pre.adM[y]

#Cell prob for Adult Males Pre-Harvest
pr.am.shot[2*y,y] <- 1 * f.pre.adM[y]
  
  #Cumulative probability of surviving to start of 2nd hunting season (surv)
  
  #Surv prob for Juv Females Post-Harvest
  surv.jf[2*y-1,y] <- s.jf.sum[y] * s.af.win[y] * s.af.sum[y+1] 

  #Surv prob for Juv Females Pre-Harvest
  surv.jf[2*y,y] <- s.jf.win[y] * s.jf.sum[y+1]

  #Surv prob for Adult Females Post-Harvest
  surv.af[2*y-1,y] <- s.af.sum[y] * s.af.win[y] * s.af.sum[y+1]

  #Surv prob for Adult Females Pre-Harvest
  surv.af[2*y,y] <- s.af.win[y] * s.af.sum[y+1]

  #Surv prob for Juv Males Post-Harvest
  surv.jm[2*y-1,y] <- s.jm.sum[y] * s.am.win[y] * s.am.sum[y+1]

  #Surv prob for Juv Males Pre-Harvest
  surv.jm[2*y,y] <- s.jm.win[y] * s.jm.sum[y+1]

  #Surv prob for Adult Males Post-Harvest
  surv.am[2*y-1,y] <- s.am.sum[y] * s.am.win[y] * s.am.sum[y+1]

  #Surv prob for Adult Males Post-Harvest
  surv.am[2*y,y] <- s.am.win[y] * s.am.sum[y+1]

} # y for first diagonal

# second y loop 
for (y in 1:yrs){   
  for (j in (y+1):yrs){
    # model prob recovery = prob survives until j, recovered in j 
    
    #Combo prob for Juv Females Post-Harvest
    pr.jf.shot[2*y-1,j] <- surv.jf[2*y-1,j-1]*f.pre.adF[j]

    #Combo prob for Juv Females Pre-Harvest
    pr.jf.shot[2*y,j] <- surv.jf[2*y,j-1]*f.pre.adF[j]

    #Combo prob for Adult Females Post-Harvest
    pr.af.shot[2*y-1,j] <- surv.af[2*y-1,j-1]*f.pre.adF[j]

    #Combo prob for Adult Females Pre-Harvest
    pr.af.shot[2*y,j] <- surv.af[2*y,j-1]*f.pre.adF[j]

    #Combo prob for Juv Males Post-Harvest
    pr.jm.shot[2*y-1,j] <- surv.jm[2*y-1,j-1]*f.pre.adM[j]

    #Combo prob for Juv Males Pre-Harvest
    pr.jm.shot[2*y,j] <- surv.jm[2*y,j-1]*f.pre.adM[j]

    #Combo prob for Adult Males Post-Harvest
    pr.am.shot[2*y-1,j] <- surv.am[2*y-1,j-1]*f.pre.adM[j]

    #Combo prob for Adult Males Post-Harvest
    pr.am.shot[2*y,j] <- surv.am[2*y,j-1]*f.pre.adM[j]


    #Update model survival to next hunting season

    #Update prob for Juv Females Post-Harvest
    surv.jf[2*y-1,j] <- surv.jf[2*y-1,j-1]*s.af.win[j-1]*s.af.sum[j] 

    #Update prob for Juv Females Pre-Harvest
    surv.jf[2*y,j]   <- surv.jf[2*y,j-1]*s.af.win[j-1]*s.af.sum[j]

    #Update prob for Adult Females Post-Harvest
    surv.af[2*y-1,j] <- surv.af[2*y-1,j-1]*s.af.win[j-1]*s.af.sum[j]

    #Update prob for Adult Females Pre-Harvest
    surv.af[2*y,j]   <- surv.af[2*y,j-1]*s.af.win[j-1]*s.af.sum[j]

    #Update prob for Juv Males Post-Harvest
    surv.jm[2*y-1,j] <- surv.jm[2*y-1,j-1]*s.am.win[j-1]*s.am.sum[j] 

    #Update prob for Juv Males Pre-Harvest
    surv.jm[2*y,j]   <- surv.jm[2*y,j-1]*s.am.win[j-1]*s.am.sum[j]

    #Update prob for Adult Males Post-Harvest
    surv.am[2*y-1,j] <- surv.am[2*y-1,j-1]*s.am.win[j-1]*s.am.sum[j]

    #Update prob for Adult Males Post-Harvest
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
bugs.data <- list(US.Harvest=US.harvest, Can.Harvest=Can.harvest, US.wings=US.wings, Can.wings=Can.wings, banded=banded, 
                  marr.juvF.shot.R=marr.juvF.shot, marr.adF.shot.R=marr.adF.shot, marr.combF.shot.R=marr.combF.shot, 
                  marr.juvM.shot.R=marr.juvM.shot, marr.adM.shot.R=marr.adM.shot, marr.combM.shot.R=marr.combM.shot, 
                  marr.ad_unkF.shot=marr.ad_unkF.shot, marr.ad_unkM.shot=marr.ad_unkM.shot,
                  marr.juvF.shot=marr.juvF.shot, marr.juvM.shot=marr.juvM.shot,
                  rel.jf=rowSums(marr.juvF.shot), rel.af=rowSums(marr.ad_unkF.shot),
                  rel.jm=rowSums(marr.juvM.shot), rel.am=rowSums(marr.ad_unkM.shot), yrs = 20,
                  US.known.age = US.wings[,7]+US.wings[,8], US.sexed.juv = US.wings[,4]+US.wings[,5],
                  US.sexed.ad = US.wings[,1]+US.wings[,2], Can.known.age = Can.wings[,7]+Can.wings[,8],
                  Can.sexed.juv = Can.wings[,4]+Can.wings[,5], Can.sexed.ad = Can.wings[,1]+Can.wings[,2])


# Initial values (not all priors listed)
inits <- function(){list(pi.US.juv.mu = runif(1, -1, 1), pi.US.juvF.mu = runif(1, -1, 1), pi.US.adF.mu = runif(1, -1, 1),
                          pi.US.juv.sd = runif(1, 0, 1), pi.US.juvF.sd = runif(1, 0, 1), pi.US.adF.sd = runif(1,0,1), 
                          pi.Can.juv.mu = runif(1, -1, 2), pi.Can.juvF.mu = runif(1, -1, 1), pi.Can.adF.mu = runif(1, -1, 1),
                          pi.Can.juv.sd = runif(1, 0, 1), pi.Can.juvF.sd = runif(1, 0, 1), pi.Can.adF.sd = runif(1,0,1),
                          s.af.sum.mu = runif(1,0,3), s.am.sum.mu = runif(1,0,3), s.jf.sum.mu = runif(1,0,3), s.jm.sum.mu = runif(1,0,3),
                          s.af.win.mu = runif(1,0,3), s.am.win.mu = runif(1,0,3), s.jf.win.mu = runif(1,0,3), s.jm.win.mu = runif(1,0,3),
                          F.jf.mu = runif(1,0.5,2), F.jm.mu = runif(1,0.5,2),   
                          N.combF.spring=c(rpois(1,958.3),rep(NA,19)),
                          N.combM.spring=c(rpois(1,1033.9),rep(NA,19)))}

# Parameters monitored   
parameters <- c("N.adF.spring", "N.adM.spring", "N.adF.fall", "N.adM.fall",
                "N.juvF.spring", "N.juvM.spring", "N.juvF.fall", "N.juvM.fall",
                "s.af.sum", "s.af.win", "s.am.sum", "s.am.win", "s.jf.sum",
                "s.jf.win", "s.jm.sum", "s.jm.win")


# MCMC settings (90 sec, converges rapidly, great mixing)
ni <- 200
nt <- 1        #10,000 posterior samples
nb <- 0
nc <- 1

ABDU.IPM.out <- jagsUI(bugs.data, inits, parameters, "Lincoln.bug", n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
