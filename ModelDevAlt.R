library(statsecol)
library(jagsUI)
library(tidyverse)
library(zoo)
data("wildebeest")

#exploratory mean-variance plot
#ggplot(testDat, aes(x = Nhat, y = sehat)) + geom_point()

#for the model:
#𝑁0 ∼ uniform(0,500k)
#𝑁𝑡|𝑁𝑡−1 ∼ Poissonl(𝜆𝑡* (𝑁𝑡−1 - c𝑡−1))
#log(𝜆𝑡) = 𝛽0 + 𝛽1R𝑡
#𝑦|𝑁 ∼ N(N𝑡, σ𝑡)

#alternative form explored in another model
#𝑦|𝑁 ∼ binomial(𝑁,𝑝) 𝑡 𝑡 𝑡

#CV is coefficient of variance: SE/E[X]
#For the censuses in the 1960s we set the CV equal to 0.3, which is about twice the average CV from the censuses in the 1970s and 1980s and reflects a lower confidence in the early censuses (Sinclair, personal communication). The standard deviations for the censuses in the 1970s were derived from asymptotic normal theory, and in keeping with this we assume that the census data are normally distributed.

#binomial variance is n*p*(1-p)
#if initial CI assumptions are correctly formulated, then these should be p and q
#interimDat <- wildebeest %>% mutate(qhat = (1 - sqrt(1 - 4*(sehat^2)/(Nhat)))/2,
#                                    phat = (1 + sqrt(1 - 4*(sehat^2)/(Nhat)))/2)

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
wildebeestAlt <- wildebeestImpute %>% mutate(Nhat = Nhat * 1000,
                                             sehat = sehat * 1000,
                                             lci = lci * 1000,
                                             uci = uci * 1000,
                                             Catch = Catch * 1000)

sink("wildebeestSSM1.txt")
cat("
model{
  #PRIORS
  n0 ~ dunif(0,0.6) #no more than 500k individuals observed before 1970
  N[1] <- n0
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01) 
  
  # Likelihood - State process
  for (t in 2:nYrs) {
    #x is rain; add a lagged term maybe?
    log(lambda[t]) <- beta0 + beta1 * X[t]
    N[t] ~ dpois(lambda[t]*(N[t-1] - c[t-1])) #subtract removals last
  }
  
  # Likelihood - Observation process
  for(tFilt in validYrs) { #why are we leaving out the last year?
    #Index out of range taking subset of  y
    y[tFilt] ~ dnorm(N[tFilt], obsTau[tFilt])
  }
  
  #DERIVED 
  #deltaN <- N[nYrs] - N[1] #total change in true population if interested
  #think of what other summary statistics we might want estimated
}", fill = TRUE)
sink()

wildebeestData <- list(nYrs = numYears,
                       validYrs = validObs[1:length(validObs)-1],
                       obsTau = 1/(wildebeestImpute$sehat^2),
                       y = wildebeestImpute$Nhat, 
                       c = wildebeestImpute$Catch,
                       X = wildebeestImpute$rain)

wildebeestInits <- function() {
  list(N = wildebeestImpute$Nhat, #abundance starting values
       beta0 = runif(1,0,1),   #given the log transform this value best to start off small
       beta1 = runif(1,0,1)) #,  #given the log transform this value best to start off small
}

wildebeestParams <- c("beta0", "beta1", "N")#, "p")

nt <- 1
nc <- 3
nb <- 10000
ni <- 100000 + nb

wildebeestFit <- jags(data = wildebeestData,
                      inits = wildebeestInits,
                      parameters.to.save = wildebeestParams,
                      model.file = "wildebeestSSM1.txt",
                      n.chains = nc,
                      n.iter = ni,
                      n.burnin = nb,
                      n.thin = nt)

for(tFilt in validObs){
  print(tFilt)
}