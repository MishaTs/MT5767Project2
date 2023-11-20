# necessary dependencies

library(tidyverse) # data wrangling

library(zoo) # back-filling

library(jagsUI) # writing BUGS language in R

library(MCMCvis) # MCMC visualization and evaluation 

library(statsecol) # load data 

library(ggplot2)

data("wildebeest")

# checking for NA values in estimate of counts, Nhat column

which(!is.na(wildebeest$Nhat), arr.ind=TRUE)

# backfilling

wildebeest <- na.locf(wildebeest, na.rm = FALSE)

# backfilling the first row 

wildebeest[1,2:6] <- wildebeest[2,2:6]

# checking the head 

head(wildebeest)

# N0 ~ uniform(0, U)
# Nt | Nt-1 ~ normal [lambda*(Nt-1 - ct-1), sigmaN]
# yt | Nt ~ normal (Nt, sigmaY)


# writing the model in BUGS
# model specification

sink("wildebeestSSMMA.txt")
cat("
model{
# PRIORS
  N1 ~ dunif(0, 0.7) # no more than 500k individuals observed before 1970
  N.est[1] <- N1
  sigmaN  ~ dunif(0, 10)             # prior for variance of state process
  tau.sigmaN  <- pow(sigmaN, 2)       
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)

# Likelihood - State process
  for (t in 1:(nyrs-1)){
  N.est[t+1] ~ dnorm(lambda[t]*(N.est[t] - c[t]), tau.sigmaN)
  lambda[t] <- beta0 + beta1 * R[t]
  }

# Likelihood - Observation process
  for(t in 1:nyrs){
  y[t] ~ dnorm(N.est[t], tau.obs[t])
  }
}", fill = TRUE)
sink()


# JAGS package data

jagsdat <- list(nyrs = nrow(wildebeest),
                       tau.obs = wildebeest$sehat^-2,
                       y = wildebeest$Nhat, 
                       c = wildebeest$Catch,
                       R = wildebeest$rain)


# set initial values for the unknown parameters

inits <- function() {
  list(
       N.est = wildebeest$Nhat, 
       sigmaN = runif(1, 0, 0.5),  
       beta0 = runif(1,-2, 2),   
       beta1 = runif(1,-2, 2))   
}

# parameters monitoring 

parms <- c("beta0", "beta1", "sigmaN", "lambda", "N.est")

# setting MCMC values 

nt <- 1 # thinning rate to reduce autocorrelation

nc <- 3 # number of chains

nb <- 2000 # burn-ins / also, called warm-ups

ni <- 70000 + nb # total iterations 


# run JAGS

out <- jags(data = jagsdat, # jags data package
            inits = inits,  # initial values
            parameters.to.save = parms, # parameters to look into 
            model.file = "wildebeestSSMMA.txt", # model file
            n.chains = nc, # number of chains
            n.iter = ni,   # number of iterations
            n.burnin = nb, # number of burn-ins
            n.thin = nt)   # thinning rate 

# trace-plots

MCMCtrace(out,                                  #the fitted model
          params = parms[1:2],                  #out parameters of interest
          type= "trace",                        #plot chain specific density curves
          iter = ni,                            #plot all iterations
          pdf = FALSE)                          #DON'T write to a PDF

MCMCtrace(out,                                  #the fitted model
          params = parms[1:2],                  #out parameters of interest
          type= "density",                      #plot chain specific density curves
          iter = ni,                            #plot all iterations
          pdf = FALSE) 


MCMCsummary(out,                 #the fitted model
            params = parms,      #out parameters of interest
            digits = 3)

#access the posterior of each parameter
mu.lambda <- out$sims.list$mu.lambda
sig.lambda <- out$sims.list$sig.lambda
sig.obs <- out$sims.list$sig.obs

Nhat <- out$sims.list$N.est
sig.lambda <- out$sims.list$lambda
dim(Nhat)

lambda_df <- data.frame(Year = 1:(nrow(wildebeest)-1),
                        Mean = out$mean$lambda,
                        Lower = out$q2.5$lambda,
                        Upper = out$q97.5$lambda)


gbars <- ggplot(data=lambda_df) +
  geom_line(aes(x=Year, y=1), color="blue", linewidth=1) +
  geom_errorbar(aes(x=Year, y=Mean, ymin=Lower, ymax = Upper, width=0)) + 
  geom_point(aes(x=Year, y=Mean), size=2.5) + 
  ylab(expression(paste("Growth rate (",lambda,")"))) + theme_bw()

gribn <- ggplot(data=lambda_df) +
  geom_line(aes(x=Year, y=1), color="blue", linewidth=1) +
  geom_ribbon(aes(x=Year, y=Mean, ymin=Lower, ymax = Upper), 
              fill="black", alpha=0.1) + 
  geom_line(aes(x=Year, y=Mean)) + 
  geom_point(aes(x=Year, y=Mean)) + 
  ylab("") + theme_bw()

cowplot::plot_grid(gbars,gribn,nrow=1)








