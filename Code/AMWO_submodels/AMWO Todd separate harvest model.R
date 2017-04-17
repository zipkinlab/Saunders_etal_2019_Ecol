#--------------------------------------------------------------------------------------------------------------
# Estimating annual harvest of AMWO during duck stamp years using annual sales as a covariate
#--------------------------------------------------------------------------------------------------------------

rm(list=ls(all=TRUE)) # clear memory

# load R2WinBUGS, set working directory for data and output, and identify WinBUGS location for your computer
library(jagsUI)
#setwd("C:/Users/arnol065/Documents/Publications/AMWO IPM/Data") 

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
  HIP[i,1] <- harvest[i,13] 
  HIP[i,2] <- harvest[(i+51),13] 
  DSS[i,1] <- harvest[i,12] 
  DSS[i,2] <- harvest[(i+51),12] 
  stamps[i,1] <- harvest[i,22]
  stamps[i,2] <- harvest[(i+51),22]}

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

# specify jags model to predict total annual harvest
sink("AMWO.harvest.bug")
cat("
    model {
  # [,1] is central pop, [,2] is eastern population
  # prior for year 1 Harvest (reverse time, [1,]=2013)
    Harvest[1,1] ~ dnorm(436095.6,1.822e-10)               # mean, 1/sd2 from 2004-2013 data
    Harvest[1,2] ~ dnorm(161459.1,9.551e-10)

  # priors for annual variation in true harvest (sigma.H) 
  # and estimated harvest (sd.HIP, sd.DSS)
  for (k in 1:2){
  sigma.H[k] ~ dunif(1,200000)
  tau.H[k] <- pow(sigma.H[k],-2)
  sd.HIP[k] ~ dunif(1,200000)         
  prec.HIP[k] <- pow(sd.HIP[k],-2)
  sd.DSS[k] ~ dunif(1,200000)
  prec.DSS[k] <- pow(sd.DSS[k],-2)
  }

  # State process (estimate true Harvest in years 2-51)
    for (k in 1:2){
      for (i in 2:yrs){
      eps[i,k] ~ dnorm(0,tau.H[k])
      Harvest[i,k] <- Harvest[i-1,k] + eps[i,k] # autoregressive model
      }
    }

# observation process, recoveries and harvest data
    for (k in 1:2){
      for (i in 1:15){
        HIP[i,k] ~ dnorm(Harvest[i,k],prec.HIP[k])}           # end HIP years
      for (j in 13:50){
        #frac[j,k] <- stamps[j,k]/max.hunters[k]
        DSS[j,k] ~ dnorm(Harvest[j,k]/cor[k],prec.DSS[k])}         # end DSS yrs. Replaced frac[j,k] with 2
    } # end k loop

  } # end jags model
    ",fill = TRUE)
sink()

# Bundle data  , add 1 to harvest, banding and recovery totals, per Seber's small-sample modification
cor <- c(4.1, 3.1)
bugs.data <- list(yrs=51, HIP=HIP, DSS=DSS, cor=cor) #, stamps=stamps)

# Initial values (not all priors listed)
#sd.HIP.inits <- c(NA, NA)
#sd.HIP.inits[1] <- runif(1,300000,400000)
#sd.HIP.inits[2] <- runif(1, 100000, 200000)

#inits <- function(){list(sd.HIP=sd.HIP.inits)}  #adjusted from 2200000 and 1000000

# Parameters monitored
parameters <- c("Harvest","sd.DSS","sd.HIP","sigma.H") #,"max.hunters")

# MCMC settings
ni <- 4000000
nt <- 200        
nb <- 2000000
nc <- 3

# call jags
AMWO.harvest.34corfactor.new <- jagsUI(bugs.data, inits=NULL, parameters, "AMWO.harvest.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)
print(AMWO.harvest.34corfactor.new,digits=1)
save(AMWO.harvest.34corfactor.new, file="AMWO Todd harvest model.Rda")

print(AMWO.harvest.jags$Rhat,digits=3)
AMWO.harvest.jags$parameters

### PLOTS########################################################################################################################################
?plot
plot(AMWO.harvest.jags$mean$Harvest[,2]~harvest[1:51,1],type="b",pch=17,col="red",ylim=c(0,1400000),xlab="Year",ylab="Harvest")
par(new=TRUE)
plot(AMWO.harvest.jags$mean$Harvest[,1]~harvest[1:51,1],type="b",pch=16,col="blue",ylim=c(0,1400000),xlab="",ylab="")
points(harvest[1:51,13]~harvest[1:51,1],col="blue")
points(harvest[1:51,12]~harvest[1:51,1],col="blue")
points(harvest[52:102,13]~harvest[1:51,1],col="red")
points(harvest[52:102,12]~harvest[1:51,1],col="red")

harvest[1,]
# plotting code (lifted from Kery & Schaub, chapter 5)
fitted <- lower <- upper <- numeric()
year <- c(1963:2013)
n.years <- 51
for (i in 1:n.years){
  fitted[i] <- mean(AMWO.harvest.jags$sims.list$Harvest[,2,i])
  lower[i] <- quantile(AMWO.harvest.jags$sims.list$Harvest[,2,i], 0.025)
  upper[i] <- quantile(AMWO.harvest.jags$sims.list$Harvest[,2,i], 0.975)}
m1 <- min(c(fitted, lower), na.rm = TRUE)
m2 <- max(c(fitted, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years),main="Eastern Harvest", ylab = "", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(fitted, type = "l", col = "blue", lwd = 2)

# plot Central Harvest
fitted <- lower <- upper <- numeric()
year <- c(1963:2013)
n.years <- 51
for (i in 1:n.years){
  fitted[i] <- mean(AMWO.harvest.jags$sims.list$Harvest[,1,i])
  lower[i] <- quantile(AMWO.harvest.jags$sims.list$Harvest[,1,i], 0.025)
  upper[i] <- quantile(AMWO.harvest.jags$sims.list$Harvest[,1,i], 0.975)}
m1 <- min(c(fitted, lower), na.rm = TRUE)
m2 <- max(c(fitted, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years),main="Eastern Harvest", ylab = "", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(fitted, type = "l", col = "blue", lwd = 2)
