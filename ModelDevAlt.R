library(statsecol)
library(jagsUI)
library(tidyverse)
data("wildebeest")

#exploratory mean-variance plot
#ggplot(testDat, aes(x = Nhat, y = sehat)) + geom_point()

#for the model:
#𝑁0 ∼ uniform(0,500k)
#𝑁𝑡|𝑁𝑡−1 ∼ Poissonl(𝜆𝑡* (𝑁𝑡−1 - c𝑡−1))
#log(𝜆𝑡) = 𝛽0 + 𝛽1R𝑡
#𝑦𝑡|𝑁𝑡 ∼binomial(𝑁𝑡,𝑝𝑡)

#CV is coefficient of variance: SE/E[X]
#For the censuses in the 1960s we set the CV equal to 0.3, which is about twice the average CV from the censuses in the 1970s and 1980s and reflects a lower confidence in the early censuses (Sinclair, personal communication). The standard deviations for the censuses in the 1970s were derived from asymptotic normal theory, and in keeping with this we assume that the census data are normally distributed.

#get row numbers for all non-na values
#manually expressed as c(2, 4, 6, 8, 12, 13, 18, 19, 23, 25, 27, 30)
validObs <- which(!is.na(wildebeest$Nhat), arr.ind=TRUE)

#binomial variance is n*p*(1-p)
#if initial CI assumptions are correctly formulated, then these should be p and q
wildeImpute <- wildebeest %>% mutate(qhat = (1 - sqrt(1 - 4*(sehat^2)/(Nhat)))/2,
                                 phat = (1 + sqrt(1 - 4*(sehat^2)/(Nhat)))/2)



wildeImpute$qhat[is.na(wildeImpute$qhat)] <- 1
wildeImpute$uci[is.na(wildeImpute$uci)] <- wildeImpute$uci[nrow(wildeImpute)]
wildeImpute[is.na(wildeImpute)] <- 0

numYears <- nrow(wildeImpute)
  
sink("wildebeestSSM2.txt")
cat("
model{
  #PRIORS
  n0 ~ dunif(0,0.5) #no more than 500k individuals observed before 1970, so likely sensible?
  N[1] <- n0
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01) 
  #beta2 ~ dnorm(0,0.01) 
  #p ~ dunif(0,1)
  
  # Likelihood - State process
  for (t in 2:nYrs) {
    #x is rain; add a lagged term maybe?
    log(lambda[t]) <- beta0 + beta1 * X[t] #+ beta2 * X[t-1]
    N[t] ~ dpois(lambda[t]*(N[t-1] - c[t-1])) #subtract removals last
  }
  
  # Likelihood - Observation process
  for(tFilt in validYrs) { 
    #y[tFilt] ~ norm(N[tFilt], obsTau)
    #binomial fit fails
    y[tFilt] ~ dbin(p[tFilt], N[tFilt])
  }
  
  #DERIVED 
  #deltaN <- N[nYrs] - N[1] #total change in true population if interested
  #think of what other summary statistics we might want estimated
}", fill = TRUE)
sink()

wildebeestData2 <- list(nYrs = numYears,
                       validYrs = validObs,
                       p = wildeImpute$phat,
                       y = wildeImpute$Nhat, 
                       c = wildeImpute$Catch,
                       X = wildeImpute$rain)

wildebeestInits2 <- function() {
  list(#p = runif(1,0,1), #probability starting value
    beta0 = runif(1,0,1),   #given the log transform this value best to start off small
    beta1 = runif(1,0,1)) #,  #given the log transform this value best to start off small
  #beta2 = runif(1,0,100))    #random number between 0 and 100
}

wildebeestParams2 <- c("beta0", "beta1", "N", "y")#, "p")

nt <- 1
nc <- 3
nb <- 10000
ni <- 100000 + nb

wildebeestFit2 <- jags(data = wildebeestData2,
                       inits = wildebeestInits2,
                       parameters.to.save = wildebeestParams2,
                       model.file = "wildebeestSSM2.txt",
                       n.chains = nc,
                       n.iter = ni,
                       n.burnin = nb,
                       n.thin = nt)