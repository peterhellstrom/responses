################################################################################
# BUGS type II functional response
################################################################################
# Simulate data
a <- 25
b <- 6
theta <- 1
n <- 50
sigma <- 1
x <- runif(n=n, min=1, max=25)
y.hat <- a*x^theta / (b^theta + x^theta)
y <- rnorm(x, mean=y.hat, sd=sigma)

curve((a*x^theta / (b^theta + x^theta)), from=0,to=30,n=1001)
points(x,y.hat,col=1)
points(x,y,col=2,pch=16)

# nls-fit
fm2 <- nls(y ~ a * x / (b + x), start=list(a=25,b=6))

summary(fm2)
coef(fm2)

plot(x,residuals(fm2))

# Call to WinBUGS
# Load required libraries
library(R2WinBUGS)
library(coda)
library(emdbook)

# Specify initial conditions for MCMC (NOT PRIORS!!!)
inits <- list(
  list(a = 20, b = 6, sigma.y=1),
  list(a = 25, b = 5, sigma.y=3),
  list(a = 27, b = 4, sigma.y=0.5)
) # end of list

# Call to WinBUGS
b.bugs <- bugs(
      data=list(x=x, y=y, n=n),
      inits,
      parameters.to.save = c("a","b","sigma.y"),
      n.thin = 1,
      n.burnin = 15000,
      model.file = "D:/data/fr typeII.bug",
      n.chains = length(inits),
      n.iter = 40000,
      #codaPkg = FALSE,
      codaPkg = TRUE,
      debug = TRUE
) # end of bugs statement

# CAUTION!!!
# Do not forget to close WinBUGS in order to do the rest!!!

model1.coda <- read.bugs(b.bugs)
est <- summary(model1.coda) # outputs a list

e <- est$statistics # Extract mean & sd's.
e[1,1] # Get the estimate for B0

# Some plotting...

plot(model1.coda)


densityplot(model1.coda)


densityplot(model1.coda)[2]

# Setting codaPkg to False

s1 <- as.mcmc.bugs(b.bugs)
gelman.diag(s1)


################################################################################
# BUGS code
# store as "fr typeII.bug"

model {
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(y.hat[i], tau.y)
    y.hat[i] <- a*x[i] / (b + x[i])
    }
  # Priors
  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)

  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif(0,100)
} # end of model
