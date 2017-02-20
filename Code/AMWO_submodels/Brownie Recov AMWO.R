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

    ## Priors and constraints, for mean phi and f, convert to logit scale
    ## I don't think these are right (Should be dnorm), bc phi[] gets converted to logit scale and not phi.mu, phi.mu can be negative, right?
    phi.adF.sum.mu ~ dunif(0, 4)        #adult female summer survival
    phi.adF.win.mu ~ dunif(0, 4)        #adult female winter survival
    f.ad.mu ~ dunif(-4, -2)           #adult recovery prob, assuming harvest probability for males & females is =
    phi.juv.sum.mu ~ dunif(0, 4)        #juv phi, summer
    phi.juv.win.mu ~ dunif(0, 4)        #juv phi, winter
    f.juv.mu ~ dunif(-4, -2)          #juv recovery prob
    phi.adM.sum.mu ~ dunif(0, 4)        #adult male phi, summer  
    phi.adM.win.mu ~ dunif(0, 4)        #adult male phi, winter 

    ## Priors for annual SD of survival and reporting rates (not very vague priors?)
    phi.adF.sum.sd ~ dunif(0, 1.5)
    phi.juv.sum.sd ~ dunif(0, 1.5)
    phi.juv.win.sd ~ dunif(0, 1.5)       
    phi.adF.win.sd ~ dunif(0, 1.5)        
    f.juv.sd ~ dunif(0, 0.5)  
    f.ad.sd ~ dunif(0, 0.5)       
    phi.adM.sum.sd ~ dunif(0, 1.5)       
    phi.adM.win.sd ~ dunif(0, 1.5)        

    # express variances as precision (1/SD^2)
    phi.adF.sum.tau <- pow(phi.adF.sum.sd,-2)
    phi.juv.sum.tau <- pow(phi.juv.sum.sd,-2)
    phi.juv.win.tau <- pow(phi.juv.win.sd,-2)
    phi.adF.win.tau <- pow(phi.adF.win.sd,-2)
    f.juv.tau <- pow(f.juv.sd,-2)
    f.ad.tau <- pow(f.ad.sd,-2)
    phi.adM.sum.tau <- pow(phi.adM.sum.sd,-2)
    phi.adM.win.tau <- pow(phi.adM.win.sd,-2)

    ## Generate season parameters for survival
    for (t in 1:yrs){
      for (i in 1:NRegion){
        for (cc in 1:NClass){
    epsilon.phi.juv.sum[t] ~ dnorm(0,phi.juv.sum.tau)
    epsilon.phi.juv.win[t] ~ dnorm(0,phi.juv.win.tau)
    epsilon.phi.adF.sum[t] ~ dnorm(0,phi.adF.sum.tau)
    epsilon.phi.adF.win[t] ~ dnorm(0,phi.adF.win.tau)
    epsilon.f.juv[t] ~ dnorm(0,f.juv.tau)
    epsilon.f.ad[t] ~ dnorm(0,f.ad.tau)
    epsilon.phi.adM.sum[t] ~ dnorm(0,phi.adM.sum.tau)
    epsilon.phi.adM.win[t] ~ dnorm(0,phi.adM.win.tau)

    logit(phi[t,1,1,i]) <- phi.juv.sum.mu + epsilon.phi.juv.sum[t]
    logit(phi[t,2,1,i]) <- phi.juv.win.mu + epsilon.phi.juv.win[t]
    logit(phi[t,1,3,i]) <- phi.adF.sum.mu + epsilon.phi.adF.sum[t]
    logit(phi[t,2,3,i]) <- phi.adF.win.mu + epsilon.phi.adF.win[t]
    logit(f[t,1,i]) <- f.juv.mu + epsilon.f.juv[t]
    logit(f[t,cc,i]) <- f.ad.mu + epsilon.f.ad[t]    #assuming recovery rates for males and females same?
    logit(phi[t,1,2,i]) <- phi.adM.sum.mu + epsilon.phi.adM.sum[t]
    logit(phi[t,2,2,i]) <- phi.adM.win.mu + epsilon.phi.adM.win[t]
    } #t
    } #i
    } #cc
    #phi.adF.sum[21] <- 1 # dummy values for easier coding--what are these for???
    #phi.adM.sum[21] <- 1

    ## Calculate number of birds released each year--by class, season, region??
    ## Why is class index always skip juveniles (see loops below as well)?
    for (t in 1:yrs){
      for (i in 1:NRegion){
        for (cc in 2:NClass){
          for (s in 1:NSeason){
      rel.juv[t,s,1,i] <- sum(marrayAMWO[t,,s,1,i])                              
      rel.ad[t,s,cc,i] <- sum(marrayAMWO[t,,s,cc,i])
    } #t
    } #i
    } #cc
    } #s
  
    ## Define multinomial likelihoods for m-arrays -- why don't we merge this loop with the above loop for quicker computing?
    for (t in 1:yrs){
      for (i in 1:NRegion){
        for (cc in 2:NClass){
          for (s in 1:NSeason){
            marrayAMWO[t,1:(yrs+1),s,1,i] ~ dmulti(pr.juv[t,,s,1,i], rel.juv[t,s,1,i])      
            marrayAMWO[t,1:(yrs+1),s,cc,i] ~ dmulti(pr.ad[t,,s,cc,i], rel.ad[t,s,cc,i])
    } #t
    } #i
    } #cc
    } #s

    ## Define cell probabilities of m-arrays
    # Main diagonal--recovery in season [t]
    for (t in 1:yrs){
      for (i in 1:NRegion){
        for (cc in 2:NClass){
      pr.juv[t,t,1,1,i] <- phi[t,1,1,i] * f[t,1,i]
      pr.juv[t,t,2,1,i] <- 1 * f[t,1,i]
      pr.ad[t,t,1,cc,i] <- phi[t,1,cc,i] * f[t,cc,i]
      pr.ad[t,t,2,cc,i] <- 1 * f[t,cc,i]

    for (k in (t+1):yrs){ #not sure how far up to move this k loop. MTF
    # probability of surviving to start of second hunting season: these are stored values to be used in subsequent diags
    surv.juv[t,t,1,1,i] <- phi[t,1,1,i] * phi[t,2,1,i] * phi[k,1,2,i]   #note that we are using 2 here in last phi but need to figure out how to do sex mixtures
    surv.juv[t,t,2,1,i] <- phi[t,2,1,i] * phi[k,1,2,i]                      #2 in second phi is placeholder for sex
    surv.ad[t,t,1,cc,i] <- phi[t,1,cc,i] * phi[t,2,cc,i] * phi[k,1,cc,i]
    surv.ad[t,t,2,cc,i] <- phi[t,2,cc,i] * phi[k,1,cc,i]

    # second and subsequent diags    
    # model prob recovery = prob survives until k, recovered in k
    pr.juv[t,k,1,1,i] <- surv.juv[t,t,1,1,i] * f[k,2,i]
    pr.juv[t,k,2,1,i] <- surv.juv[t,t,2,1,i] * f[k,2,i]
    pr.ad[t,k,1,cc,i] <- surv.ad[t,t,1,cc,i] * f[k,cc,i]
    pr.ad[t,k,2,cc,i] <- surv.ad[t,t,2,cc,i] * f[k,cc,i]

    # survival to next hunting season (all adult survival now)
    surv.juv[t,k,1,1,i] <- surv.juv[t,t,1,1,i] * phi[t,2,2,i] * phi[k,1,2,i] #check to make sure it is not k and k + 1 MTF
    surv.juv[t,k,2,1,i] <- surv.juv[t,t,2,1,i] * phi[t,2,2,i] * phi[k,1,2,i]
    surv.ad[t,k,1,cc,i] <- surv.ad[t,t,1,cc,i] * phi[t,2,cc,i] * phi[k,1,cc,i]
    surv.ad[t,k,2,cc,i] <- surv.ad[t,t,2,cc,i] * phi[t,2,cc,i] * phi[k,1,cc,i]
    } #k

    # Left of main diag
    for (l in 1:(t-1)){
      pr.juv[t,l,1,1,i] <- 0
      pr.juv[t,l,2,1,i] <- 0
      pr.ad[t,l,1,cc,i] <- 0
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



