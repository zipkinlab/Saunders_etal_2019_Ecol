##########################################################################
# Based off of Todd's ABDU Brownie model
# 2-season band recovery model for AMWO (summer and fall/winter)
# Female and Male, Brownie
# Juvs transition to adults after first winter, no unk adults right?
# Juvs not sexed so S rates above first diagnol are mix of S-AdF & S-AdM??
##########################################################################

sink("AMWO.Brownie")
cat("
    model {

    ## Priors and constraints, for mean S and f, convert to logit scale
    s.adF.sum.mu ~ dunif(0, 4)        #adult female summer survival
    s.adF.win.mu ~ dunif(0, 4)        #adult female winter survival
    f.ad.mu ~ dunif(-4, -2)           #adult recovery prob
    s.juv.sum.mu ~ dunif(0, 4)        #juv S, summer
    s.juv.win.mu ~ dunif(0, 4)        #juv S, winter
    f.juv.mu ~ dunif(-4, -2)          #juv recovery prob
    s.adM.sum.mu ~ dunif(0, 4)        #adult male S, summer  
    s.adM.win.mu ~ dunif(0, 4)        #adult male S, winter 

    ## Priors for annual SD of survival and reporting rates
    s.adF.sum.sd ~ dunif(0, 1.5)
    s.juv.sum.sd ~ dunif(0, 1.5)
    s.juv.win.sd ~ dunif(0, 1.5)       
    s.adF.win.sd ~ dunif(0, 1.5)        
    f.juv.sd ~ dunif(0, 0.5)  
    f.ad.sd ~ dunif(0, 0.5)       
    s.adM.sum.sd ~ dunif(0, 1.5)       
    s.adM.win.sd ~ dunif(0, 1.5)        

    # express variances as precision (1/SD^2)
    s.adF.sum.tau <- pow(s.adF.sum.sd,-2)
    s.juv.sum.tau <- pow(s.juv.sum.sd,-2)
    s.juv.win.tau <- pow(s.juv.win.sd,-2)
    s.adF.win.tau <- pow(s.adF.win.sd,-2)
    f.juv.tau <- pow(f.juv.sd,-2)
    f.ad.tau <- pow(f.ad.sd,-2)
    s.adM.sum.tau <- pow(s.adM.sum.sd,-2)
    s.adM.win.tau <- pow(s.adM.win.sd,-2)

    ## Generate season parameters for survival
    for (t in 1:yrs){
      for (i in 1:NRegion){
        for (cc in 1:NClass){
    epsilon.s.juv.sum[t] ~ dnorm(0,s.juv.sum.tau)
    epsilon.s.juv.win[t] ~ dnorm(0,s.juv.win.tau)
    epsilon.s.adF.sum[t] ~ dnorm(0,s.adF.sum.tau)
    epsilon.s.adF.win[t] ~ dnorm(0,s.adF.win.tau)
    epsilon.f.juv[t] ~ dnorm(0,f.juv.tau)
    epsilon.f.ad[t] ~ dnorm(0,f.ad.tau)
    epsilon.s.adM.sum[t] ~ dnorm(0,s.adM.sum.tau)
    epsilon.s.adM.win[t] ~ dnorm(0,s.adM.win.tau)

    logit(phi[t,1,1,i]) <- s.juv.sum.mu + epsilon.s.juv.sum[t]
    logit(phi[t,2,1,i]) <- s.juv.win.mu + epsilon.s.juv.win[t]
    logit(phi[t,1,3,i]) <- s.adF.sum.mu + epsilon.s.adF.sum[t]
    logit(phi[t,2,3,i]) <- s.adF.win.mu + epsilon.s.adF.win[t]
    logit(f[t,1,i]) <- f.juv.mu + epsilon.f.juv[t]
    logit(f[t,cc,i]) <- f.ad.mu + epsilon.f.ad[t]    #assuming recovery rates for males and females same?
    logit(phi[t,1,2,i]) <- s.adM.sum.mu + epsilon.s.adM.sum[t]
    logit(phi[t,2,2,i]) <- s.adM.win.mu + epsilon.s.adM.win[t]
      } #t
    }  #i
    } #cc
    #s.adF.sum[21] <- 1 # dummy values for easier coding--what are these for???
    #s.adM.sum[21] <- 1

    ## Calculate number of birds released each year--by class, season, region??
    for (t in 1:yrs){
      for (i in 1:NRegion){
        for (cc in 2:NClass){
          for (s in 1:NSeason){
      rel.juv[t,s,1,i] <- sum(marrayAMWO[t,,s,1,i])                              
      rel.ad[t,s,cc,i]  <-sum(marrayAMWO[t,,s,cc,i])
    } #t
    } #i
    } #cc
    } #s
  
    ## Define multinomial likelihoods for m-arrays
    for (t in 1:(2*yrs)){        #2*??
      for (i in 1:NRegion){
        for (cc in 2:NClass){
          for (s in 1:NSeason){
            marrayAMWO[t,1:(yrs+1),s,1,i] ~ dmulti(pr.juv[t,,s,1,i], rel.juv[t,s,1,i])      
            marrayAMWO[t,1:(yrs+1),s,cc,i] ~ dmulti(pr.ad[t,,s,cc,i], rel.ad[t,s,cc,i])
    }}}}

    ## Define cell probabilities of m-arrays
    # Main diagonal--recovery in season [t]
    for (t in 1:yrs){
      for (i in 1:NRegion){
        for (cc in 2:NClass){
      pr.juv[t-1,t,1,1,i] <- phi[t,1,1,i] * f[t,1,i]
      pr.juv[t,t,2,1,i] <- 1 * f[t,1,i]
      pr.ad[t-1,t,1,cc,i] <- phi[t,1,cc,i] * f[t,cc,i]
      pr.ad[t,t,2,cc,i] <- 1 * f[t,cc,i]

    # probability of surviving to start of second hunting season: these are stored values to be used in subsequent diags
    surv.juv[t-1,t,1,1,i] <- phi[t-1,1,1,i] * phi[t-1,2,1,i] * phi[t,1,2,i]   #note that we are using 2 here in last phi but need to figure out how to do sex mixtures
    surv.juv[t,t,2,1,i] <- phi[t-1,2,1,i] * phi[t,1,2,i]                      #2 in second phi is placeholder for sex
    surv.ad[t-1,t,1,cc,i] <- phi[t-1,1,cc,i] * phi[t-1,2,cc,i] * phi[t,1,cc,i]
    surv.ad[t,t,2,cc,i] <- phi[t-1,2,cc,i] * phi[t,1,cc,i]

    # second and subsequent diags
    for (k in (y+1):yrs){
    # model prob recovery = prob survives until k, recovered in k
    pr.juv[t-1,k,1,1,i] <- surv.juv[t-1,k-1,1,1,i] * f[t,2,i]
    pr.juv[t,k,2,1,i] <- surv.juv[t,k-1,2,1,i] * f[t,2,i]
    pr.ad[t-1,k,1,cc,i] <- surv.ad[t-1,k-1,1,cc,i] * f[t,cc,i]
    pr.ad[t,k,2,cc,i] <- surv.ad[t,k-1,2,cc,i] * f[t,cc,i]

    # survival to next hunting season (all adult survival now)
    surv.juv[t-1,k,1,1,i] <- surv.juv[t-1,k-1,1,1,i] * phi[t-1,2,2,i] * phi[t,1,2,i]
    surv.juv[t-1,k,2,1,i] <- surv.juv[t,k-1,2,1,i] * phi[t-1,2,2,i] * phi[t,1,2,i]
    surv.ad[t-1,k,1,cc,i] <- surv.ad[t-1,k-1,1,cc,i] * phi[t-1,2,cc,i] * phi[t,1,cc,i]
    surv.ad[t-1,k,2,cc,i] <- surv.ad[t,k-1,2,cc,i] * phi[t-1,2,cc,i] * phi[t,1,cc,i]
    } #k

    # Left of main diag
    for (l in 1:(t-1)){
      pr.juv[t-1,l,1,1,i] <- 0
      pr.juv[t,l,2,1,i] <- 0
      pr.ad[t-1,l,1,cc,i] <- 0
      pr.ad[t,l,2,cc,i] <- 0
    } #l
  } #t

    # Last column: probability of non-recovery
    for (w in 1:(2*yrs)){
    pr.juv[t,(yrs+1)] <- 1-sum(pr.juv[t,1:yrs,,1,])
    pr.ad[t,(yrs+1)] <- 1-sum(pr.ad[t,1:yrs,,,])
    } #w
  } #cc
} #i
      
    } # end bugs model
    ",fill = TRUE)
sink()



