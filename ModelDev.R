library(statsecol)
library(jagsUI)
library(tidyverse)
data("wildebeest")

#for the model:
#𝑁0 ∼uniform(0,𝑈)
#𝑁𝑡|𝑁𝑡−1 ∼ Poissonl(𝜆𝑡𝑁𝑡−1)
#log(𝜆𝑡) = 𝛽0 + 𝛽1𝑋𝑡
#𝑦|𝑁 ∼binomial(𝑁,𝑝) 𝑡 𝑡 𝑡

#CV is coefficient of variance: SE/E[X]
#For the censuses in the 1960s we set the CV equal to 0.3, which is about twice the average CV from the censuses in the 1970s and 1980s and reflects a lower confidence in the early censuses (Sinclair, personal communication). The standard deviations for the censuses in the 1970s were derived from asymptotic normal theory, and in keeping with this we assume that the census data are normally distributed.


#binomial variance is n*p*(1-p)
#if initial CI assumptions are correctly formulated, then these should be p and q
testDat <- wildebeest %>% mutate(qhat = (1 - sqrt(1 - 4*(sehat^2)/(Nhat)))/2,
                                 phat2 = (1 + sqrt(1 - 4*(sehat^2)/(Nhat)))/2)

#get row numbers for all non-na values
validObs <- which(!is.na(wildebeest$Nhat), arr.ind=TRUE)

ggplot(testDat, aes(x = Nhat, y = sehat)) + geom_point()


numYears <- nrow(wildebeest)
N <- Y <- numeric(years) 
X <- runif(years-1,-2,2) #scaled covariate 
beta_0 <- log(1.02) #intercept 
beta_1 <-0.05 #slope
lambda <-exp(beta_0 + beta_1*X) #lambda ~ X 
p <- 0.7 #false neg 
N[1]<- 50 #initial N 
Y[1]<- rbinom(1,N[1],p) #initial Y

sink("wildebeestSSM1.txt")
cat("
model{
  #PRIORS
  n0 ~ dunif(0,0.5)
  N.est[1] <- n0
  beta0 ~ dnorm(0,0.01)
  beta1~dnorm(0,0.01) 
  p ~ dunif(0,1)
  
  # Likelihood - State process
  for (t in 2:(nyrs-1)){
    log(lambda[t-1]) <- beta0 + beta1 * X[t-1] #x is rain
    N.est[t] ~ dpois(lambda[t-1]*N.est[t-1])
  }
  
  # Likelihood - Observation process
  for(t in 1:validYrs){
    y[i] ~ dbin(p, N.est[t])
  }
  
  #DERIVED
  deltaN <- N[nyrs] - N[1]
}", fill = TRUE)
sink()

wildedata <- list(nyrs = numYears,
                  validYrs = validObs,
                  y = wildbeest$Nhat, 
                  c = wildebeest$Catch,
                  X = wildebeest$rain)