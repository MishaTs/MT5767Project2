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
#data{
#  logy <- log(y)
#}
model{

  # Priors and constraints
  log.n1 ~ dnorm(-1.3, 0.01)
  logN.est[1] <- log.n1 
  #mu.r ~ dnorm(1, 0.001) 
  sig.r ~ dunif(0, 1) 
  sig2.r <- pow(sig.r, 2)
  tau.r <- pow(sig.r, -2)
  sig.obs ~ dunif(0, 1) # Prior for sd of obs. process
  sig2.obs <- pow(sig.obs, 2)
  tau.obs <- pow(sig.obs, -2)
  beta0 ~ dnorm(0,0.001) 
  beta1 ~ dnorm(0,0.001)

  # Likelihood - State process
  for(t in 1:(nyrs-1)){
    mu.r[t] <- beta0 + beta1*R[t]
    r[t] ~ dnorm(mu.r[t], tau.r)
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
wildedata <- list(y = log(wildebeest$Nhat), 
                  nyrs = nrow(wildebeest), 
                  validYrs = validObs,
                  R = wildebeest$rain)

# Initial values
wildeinits <- function(){
  list(sig.r = runif(1, 0, 1), 
       beta0 = rnorm(1),
       beta1 = rnorm(1),
       #mu.r = rnorm(1),
       sig.obs = runif(1, 0, 1))
}

# Parameters monitored
wildeparms <- c("beta0", "beta1", "sig.obs", "sig.r", "r", "N.est")

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
          params = wildeparms[1:4], #out parameters of interest
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          ind = TRUE)                #chain specific densities

MCMCsummary(wildeout,
            params = wildeparms[1:4]) #out parameters of interest

wilde_traj <- data.frame(Year = wildebeest$year,
                         Mean = wildeout$mean$N.est,
                         Lower = wildeout$q2.5$N.est,
                         Upper = wildeout$q97.5$N.est,
                         Obs = wildebeest$Nhat,
                         LowerObs = wildebeest$lci,
                         UpperObs = wildebeest$uci)


ggplot(data = wilde_traj) + 
  geom_ribbon(aes(x=Year, y=Mean, ymin=Lower, ymax=Upper),
              fill="cyan", alpha = 0.25) +
  geom_line(aes(x=Year, y=Mean), linewidth=1, color="blue") + 
  geom_point(aes(x=Year, y=Obs), size=1.2) +
  geom_line(data = na.omit(wilde_traj), aes(x=Year, y=Obs)) +
  geom_errorbar(aes(x=Year, 
                    y=Obs,
                    ymin=LowerObs,
                    ymax=UpperObs), width=0, color="grey") +
  theme_bw()