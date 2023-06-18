# BUGS type II functional response
# 2 prey types

# works... but should be rewritten.
# How to incoporate a grouping variable?
# but rewrite to JAGS code.

# Simulate some data
a <- c(25,50)
b <- c(6,9)
n <- c(50,30)
sigma <- c(1,2)

x <- lapply(1:2, function(i) runif(n=n[i], min=1, max=25))
y.hat <- lapply(1:2, function(i) a[i]*x[[i]] / (b[i]+x[[i]]))
y <- lapply(1:2, function(i) rnorm(x[[i]],mean=y.hat[[i]],sd=sigma[i]))

n <- list(
n[1],
n[2]
)

dev.new()
plot(x=c(0,30),y=c(0,max(unlist(y))),type="n",xlab="Prey density",ylab="Prey eaten",las=1,font.lab=2)
for (i in 1:2) curve((a[i]*x / (b[i]+x)), from=0,to=30,n=1001, add=T)
for (i in 1:2) points(x[[i]],y.hat[[i]])
for (i in 1:2) points(x[[i]],y[[i]],col=2,pch=16)

# nls-fit
for (i in 1:2) {
	y.new <- y[[i]]
	x.new <- x[[i]]
	fm1 <- nls(y.new~a*x.new / (b+x.new), start=list(a=25,b=6))
	print(summary(fm1))
	print(coef(fm1))
}

x

x1 <- x[[1]]
x2 <- x[[2]]
y1 <- y[[1]]
y2 <- y[[2]]
n1 <- n[[1]]
n2 <- n[[2]]

################################################################################
# Call to WinBUGS
# Load required libraries
library(R2WinBUGS)
library(coda)
library(emdbook)

# Specify initial conditions for MCMC (NOT PRIORS!!!)
inits <- list(
  list(a1 = 25, a2 = 40, b1 = 10, b2 = 12, sigma.y1 = 1, sigma.y2 = 1),
  list(a1 = 20, a2 = 60, b1 = 12, b2 = 8, sigma.y1 = 2, sigma.y2 = 2),
  list(a1 = 10, a2 = 30, b1 = 10, b2 = 10, sigma.y1 = 3, sigma.y2 = 3)
  ) # end of list

# Call to WinBUGS
b.bugs <- bugs(
      data=list(x1=x1, x2=x2, y1=y1, y2=y2, n1=n1, n2=n2),
      inits,
      parameters.to.save = c("a1","a2","b1","b2","sigma.y1","sigma.y2","diff.a","diff.b","diff.sigma.y"),
      n.thin = 1,
      n.burnin = 15000,
      model.file = "D:/data/fr typeII two prey.bug",
      n.chains = length(inits),
      n.iter = 40000,
      #codaPkg = FALSE,
      codaPkg = TRUE,
      debug = TRUE
) # end of bugs statement

# Do not forget to close WinBUGS in order to do the rest!!!

model1.coda <- read.bugs(b.bugs)
est <- summary(model1.coda) # outputs a list

e <- est$statistics # Extract mean & sd's.
e[1,1] # Get the estimate for B0

e

# Some plotting...

plot(model1.coda)


densityplot(model1.coda)


densityplot(model1.coda)[6]

# Setting codaPkg to False

s1 <- as.mcmc.bugs(b.bugs)
gelman.diag(s1)


################################################################################
# BUGS code
# store as "fr typeII two prey.bug"

model {
  # Likelihood
  
  for (i in 1:n1) {
    y1[i] ~ dnorm(y.hat1[i], tau.y1)
    y.hat1[i] <- a1*x1[i] / (b1 + x1[i])
    }
  
  for (i in 1:n2) {
    y2[i] ~ dnorm(y.hat2[i], tau.y2)
    y.hat2[i] <- a2*x2[i] / (b2 + x2[i])
    }

  # Priors
  a1 ~ dnorm(0, 0.0001)
  b1 ~ dnorm(0, 0.0001)
  
  a2 ~ dnorm(0, 0.0001)
  b2 ~ dnorm(0, 0.0001)
  
  tau.y1 <- pow(sigma.y1, -2)
  sigma.y1 ~ dunif(0,100)
  
  tau.y2 <- pow(sigma.y2, -2)
  sigma.y2 ~ dunif(0,100)

# Calculated variables
diff.a <- a1 - a2
diff.b <- b1 - b2
diff.sigma.y <- sigma.y1 - sigma.y2

} # end of model
