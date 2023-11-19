library(statsecol)
library(jagsUI)
library(tidyverse)
library(zoo)
library(MCMCvis)
data("wildebeest")

#wildeTrans <- wildebeest %>% mutate(logy = log(Nhat),
#                                    logse = sqrt(log(sehat^2/Nhat^2 + 1)))

#Specify model in BUGS language
sink("wildessm.txt")
cat("
model{

  # Priors and constraints
  log.n1 ~ dnorm(-1.3, 0.01)
  logN.est[1] <- log.n1 
  mu.r ~ dnorm(1, 0.001) 
  sig.r ~ dunif(0, 1) 
  sig2.r <- pow(sig.r, 2)
  tau.r <- pow(sig.r, -2)
  sig.obs ~ dunif(0, 1) # Prior for sd of obs. process
  sig2.obs <- pow(sig.obs, 2)
  tau.obs <- pow(sig.obs, -2)

  # Likelihood - State process
  for(t in 1:(nyrs-1)){
    r[t] ~ dnorm(mu.r, tau.r)
    logN.est[t+1] <- logN.est[t] + r[t]
  }

  # Likelihood - Observation process
  for (t in validYrs) {
    y[t] ~ dnorm(log(N.est[t]), tau.obs)
  }
  
  # Derive population sizes on real scale
  for (t in 1:nyrs) {
    N.est[t] <- exp(logN.est[t])
  }
}
",fill = TRUE)
sink()

# Bundle data
wildedata <- list(y = log(wildebeestImpute$Nhat), 
                   nyrs = nrow(wildebeestImpute), 
                   validYrs = validObs)

# Initial values
wildeinits <- function(){
  list(sig.r = runif(1, 0, 1), 
       mu.r = rnorm(1),
       sig.obs = runif(1, 0, 1))
}

# Parameters monitored
wildeparms <- c("mu.r", "sig.obs", "sig.r", "r", "N.est")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3

wildeout <- jags(data = wildedata,
                  inits = wildeinits,
                  parameters.to.save = wildeparms,
                  model.file = "wildessm.txt",
                  n.chains = nc,
                  n.iter = ni,
                  n.burnin = nb,
                  n.thin = nt)

MCMCtrace(wildeout,                 #the fitted model
          params = wildeparms[1:3], #out parameters of interest
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind = TRUE)                #chain specific densities

MCMCsummary(wildeout,
            params = wildeparms[1:3]) #out parameters of interest