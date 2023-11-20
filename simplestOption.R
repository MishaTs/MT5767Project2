library(statsecol)
library(jagsUI)
library(tidyverse)
library(zoo)
library(MCMCvis)
data("wildebeest")

#get row numbers for all non-na values
#manually expressed as c(2, 4, 6, 8, 12, 13, 18, 19, 23, 25, 27, 30)
validObs <- which(!is.na(wildebeest$Nhat), arr.ind=TRUE)
numYears <- nrow(wildebeest)

# Error in node N[19]: Invalid parent values
# imputing data solves partially
wildebeestImpute <- na.locf(wildebeest, na.rm = FALSE)
#fill in 1960 with 1961 values too, since there's no better approximation
#manually done to avoid wiping rain, catch, and year values
wildebeestImpute[1,2:6] <- wildebeestImpute[2,2:6]

#Specify model in BUGS language
sink("wildessmBasic.txt")
cat("
model{

  # Priors and constraints
  #very clean at 0.6, but our y > N when set to 1. Set to a happy middle ground of 0.7 but I'm flexible
  n1 ~ dunif(0,0.7) #the model is HIGHLY sensitive to this prior
  N.est[1] <- round(n1)
  beta0 ~ dnorm(0,0.001) 
  beta1 ~ dnorm(0,0.001)
  sig.r ~ dunif(0, 1) 
  sig2.r <- pow(sig.r, 2)
  tau.r <- pow(sig.r, -2)

  # Likelihood - State process
  for(t in 1:(nyrs-1)){
    log.lambda[t] <- beta0 + beta1*R[t]
    log(lambda[t]) <- log.lambda[t]
    N.est[t+1] ~ dnorm(lambda[t]*(N.est[t] - c[t]),tau.r)
  }

  # Likelihood - Observation process
  for (t in validYrs) {
    y[t] ~ dnorm(N.est[t], obs.tau[t])
  }
  
}
",fill = TRUE)
sink()

# Bundle data
# Works with both imputed and raw data
wildedata <- list(y = wildebeest$Nhat, 
                  nyrs = nrow(wildebeest), 
                  validYrs = validObs,
                  R = wildebeest$rain,
                  obs.tau = wildebeest$sehat^-2,
                  c = wildebeest$Catch)

# Initial values
wildeinits <- function(){
  list(beta0 = rnorm(1),
       beta1 = rnorm(1),
       sig.r = runif(1))#,
       #N.est = wildebeestImpute$Nhat*1000)
}

# Parameters monitored
wildeparms <- c("beta0", "beta1", "sig.r", "lambda", "N.est")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3

wildeout <- jags(data = wildedata,
                 inits = wildeinits,
                 parameters.to.save = wildeparms,
                 model.file = "wildessmBasic.txt",
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




#project values forward 5 years from 1990-1994
nproj <- 5

#assume that illegal harvesting continues at current levels
#use average observed rainfall for future projections
#impute last year values for Nhat and sehat; they aren't referenced though
wildedata_proj <- list(y = c(wildebeestImpute$Nhat, rep(wildebeestImpute$Nhat[nrow(wildebeestImpute)], nproj)), 
                       nyrs = nrow(wildebeestImpute) + nproj, 
                       validYrs = validObs,
                       R = c(wildebeestImpute$rain, rep(mean(wildebeestImpute$rain), nproj)),
                       obs.tau = c(wildebeestImpute$sehat^-2, rep(wildebeestImpute$sehat[nrow(wildebeestImpute)]^-2, nproj)),
                       c = c(wildebeestImpute$Catch, rep(wildebeestImpute$Catch[nrow(wildebeestImpute)], nproj)))

wildeproj <- jags(data = wildedata_proj,
                  inits = wildeinits,
                  parameters.to.save = wildeparms,
                  model.file = "wildessmBasic.txt",
                  n.chains = nc,
                  n.iter = ni,
                  n.burnin = nb,
                  n.thin = nt)

wilde_proj <- data.frame(Year = c(wildebeest$year, 1990:1994),
                         Mean = wildeproj$mean$N.est,
                         Lower = wildeproj$q2.5$N.est,
                         Upper = wildeproj$q97.5$N.est,
                         Obs = c(wildebeest$Nhat,rep(NA,nproj)),
                         LowerObs = c(wildebeest$lci,rep(NA,nproj)),
                         UpperObs = c(wildebeest$uci,rep(NA,nproj)))

ggplot(data = wilde_proj) + 
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