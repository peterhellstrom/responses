##################################################
# EXPONENTIAL DISTRIBUTION
##################################################

##################################################
# Plot probability density function
dev.new()
plot(x=c(0,3), y=c(0,1), xlab="", ylab="", type="n")
curve(dexp(x=x, rate=1), col="blue", add=T)
curve(dexp(x=x, rate=1/2), col="green", add=T)
curve(dexp(x=x, rate=1/4), col="red", add=T)

# Plot cumulative density function
dev.new()
plot(x=c(0,3), y=c(0,1), xlab="", ylab="", type="n")
curve(pexp(q=x, rate=1), col="blue", add=T)
curve(pexp(q=x, rate=1/2), col="green", add=T)
curve(pexp(q=x, rate=1/4), col="red", add=T)

# Probability density function
lambda <- 1 # rate

xv <- seq(0,3,0.1)
lambda*exp(-lambda*xv) # density defined by this expression
dexp(x=xv, rate=lambda) # R-function, should give exactly the same result

# The mean is equal to 1/rate.

n <- 100000
lambda <- 0.5
.x <- rexp(n=n, rate=lambda)
hist(.x, breaks=30, col="lightgrey", freq=FALSE, main="Sample from exponential distribution")

mean(.x); 1/lambda

##################################################
# Generate random deviates
n <- 50000
lambda <- 1/4

.x <- rexp(n=n, rate=lambda)
range(.x)

# Create a histogram and draw the probability density function over the bars.
dev.new()
hist(.x, breaks=30, col="lightgrey", freq=FALSE, main="Sample from exponential distribution")
curve(dexp(x=x, rate=lambda), col="blue", lwd=2, add=T)

##################################################
# Estimate moments of a sample from the exponential distribution
# Requires 
library(fitdistrplus)
##################################################

# Create random sample from lognormal distribution
n <- 10000
lambda <- 0.25

.x <- rexp(n=n, rate=lambda)

plotdist(.x, breaks=30, col="red", type="l")

descdist(.x, boot=999)

f1.exp <- fitdist(.x, "exp") # Returns moments on log-scale
mean(.x); sd(.x) # For comparison
1/mean(.x)

summary(f1.exp)
gofstat(f1.exp)
plot(f1.exp)

# Compare with lognormal distribution
f1.lnorm <- fitdist(.x, "lnorm")
summary(f1.lnorm)
plot(f1.lnorm)

# Bootstrap
b1 <- bootdist(f1.exp, niter=100)
summary(b1)

plot(b1)

##################################################
