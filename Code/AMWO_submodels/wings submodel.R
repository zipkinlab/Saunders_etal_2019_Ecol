sink("WINGS.jags")
cat("
    model {
for(p in 1:2){
for (t in 1:yrs){
  
  pi.age[t,1,p] ~ dunif(0,1)
  #pi.sex[t,1,p] ~ dunif(0,1)
  #pi.sex[t,2,p] ~ dunif(0,1)
  
  wings.age[t,1,p] ~ dbin(pi.age[t,1,p], wings[t,p])          # pi.age is the proportion of juveniles to adults    
  #wings.sex[t,,p] ~ dmulti(pi.sex[t,,p], wings.age[t,2,p])    # this needs to stay dmulti because of pi.sex on left-side of equations below; can't put (1-pi.sex[t,1,p]) as a derived quantity    
  
  
  H[t,1,p] <- pi.age[t,1,p] * H.total[t,p]
  #pi.combo.m[t,p] <- (1-pi.age[t,1,p])*pi.sex[t,1,p]
  #H[t,2,p] ~ dbin(pi.combo.m[t,p], H.total[t,p])
  #pi.combo.f[t,p] <- (1-pi.age[t,1,p])*pi.sex[t,2,p]
  #H[t,3,p] ~ dbin(pi.combo.f[t,p], H.total[t,p])
}#t
}#p
    } # end bugs model
    ",fill = TRUE)
sink()

bugs.data <- list(yrs=dim(marrayAMWO)[1], wings.age=wings.age, wings.sex=wings.sex, H.total=H.total, wings=wings)  

wings.sex.inits <- wings.sex
wings.sex.inits[is.na(wings.sex.inits)] <- 0
wings.inits <- wings
wings.inits[is.na(wings.inits)] <- 0
wings.age.inits <- wings.age
wings.age.inits[is.na(wings.age.inits)] <- 0
H.total.inits <- H.total
H.total.inits[is.na(H.total.inits)] <- 0
inits <- function(){list(wings.sex = wings.sex.inits, H.total = H.total.inits, wings = wings.inits)}

# Parameters monitored
parameters <- c("H")  # add to this before running!

# MCMC settings
ni <- 1000
nt <- 1
nb <- 100
nc <- 1

# Run JAGS
AMWO.WINGS.jags <- jagsUI(bugs.data, inits=inits, parameters, "WINGS.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE) #changed from inits=NULL

