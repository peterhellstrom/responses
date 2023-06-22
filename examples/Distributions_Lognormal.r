##################################################
# LOG-NORMAL DISTRIBUTION
##################################################

##################################################
# Plot probability density function
plot(x=c(0,3), y=c(0,2), xlab="", ylab="", type="n")
curve(dlnorm(x=x, meanlog=0, sdlog=1), col="blue", add=T)
curve(dlnorm(x=x, meanlog=0, sdlog=1/2), col="green", add=T)
curve(dlnorm(x=x, meanlog=0, sdlog=1/4), col="red", add=T)

# Plot cumulative density function
plot(x=c(0,3), y=c(0,1), xlab="", ylab="", type="n")
curve(plnorm(q=x, meanlog=0, sdlog=1), col="blue", add=T)
curve(plnorm(q=x, meanlog=0, sdlog=1/2), col="green", add=T)
curve(plnorm(q=x, meanlog=0, sdlog=1/4), col="red", add=T)

##################################################
# Generate a random sample from a lognormal distribution.
meanlog <- 1
sdlog <- 1/2

# Note that for a log-normal distribution, the median (of a random sample x) equals exp(meanlog)
exp(meanlog)
# The mean equals exp(meanlog + (sdlog^2)/2)
exp(meanlog + (sdlog^2)/2)

##################################################
# In Excel's PopTools, you have to recalculate the input values to create a sample from a lognormal distribution (if you want to match it to R's formulation!):
# PopTools (and other functions such as LNORMINV, has shape & scale parameters on back-transformed scale as inputs!

# Functions to do the correct calculations:
# Source: Evans et al 1993 Statistical Distributions (2nd ed), p. 102-105

lnorm2norm <- function(mu,sigma) {
	names(mu) <- names(sigma) <- NULL
	 m <- exp(mu + (sigma^2)/2)
	 v <- exp(2*mu + sigma^2)*(exp(sigma^2)-1)
	 sd <- sqrt(v)
	 c(m=m, sd=sd, v=v)
}

norm2lnorm <- function(m,sd) {
	names(m) <- names(sd) <- NULL
	v <- sd^2
	mu <- log((m^2) / sqrt(v + m^2))
	sigma <- sqrt(log(v/(m^2) + 1))
	c(mu=mu,sigma=sigma,v=v)
}

##################################################
# EXAMPLES
# Create a sample with mean=1 and variance=2
pars <- norm2lnorm(1,sqrt(2))
n <- 1000000
.x <- rlnorm(n, meanlog=pars["mu"], sdlog=pars["sigma"])
mean(.x); median(.x); sd(.x); var(.x)

# Create a sample with meanlog=1 and sdlog=1/2
pars <- lnorm2norm(1,1/2)
pars
pars <- norm2lnorm(pars["m"], pars["sd"])
n <- 1000000
.x <- rlnorm(n, meanlog=pars["mu"], sdlog=pars["sigma"])
.x <- log(.x)
mean(.x); median(.x); sd(.x); var(.x)

# "Excel-style"
# Generate a sample from a normal distribution with mean=meanlog and sd=sdlog, then take exponent of the random deviates:

n <- 100000
meanlog <- 1
sdlog <- 1/2

par(mfrow=c(1,2))
.x <- rnorm(n, mean=meanlog, sd=sdlog)
range(.x)
	hist(.x, col="lightgrey", breaks=30, freq=FALSE)
	curve(dnorm(x, mean=meanlog, sd=sdlog), add=TRUE, col="blue", lwd=2)
.x <- exp(.x)
	hist(.x, col="lightgrey", breaks=30, freq=FALSE)
	curve(dlnorm(x, meanlog=meanlog, sdlog=sdlog), add=TRUE, col="blue", lwd=2)
par(mfrow=c(1,1))

mean(.x); median(.x); sd(.x); var(.x)

# Generate random deviates
n <- 50000
meanlog <- 1
sdlog <- 1/2

.x <- rlnorm(n=n, meanlog=meanlog, sdlog=sdlog)

# The mean m and standard deviation of this sample should approximate:
lnorm2norm(meanlog, sdlog)

# Check that this is roughly correct:
mean(.x); mean(log(.x))
sd(.x); sd(log(.x))

# Another way of generating random deviates:
# First draw random normal deviates with mean=0 and sd=1.
.x2 <- exp(meanlog + sdlog*rnorm(n=n, mean=0, sd=1))

# Plot random sample, together with pdf. Also plot the random sample on log-scale.
# Create a histogram and draw the probability density function over the bars.
# Then add the log-transformed values and the corresponding pdf.

par(mfrow=c(1,2))
	hist(.x, breaks=30, col="lightgrey", freq=FALSE)
	curve(dlnorm(x=x, meanlog=meanlog, sdlog=sdlog), col="blue", lwd=2, add=T)

	hist(log(.x), breaks=30, col="lightgrey", freq=FALSE)
	curve(dnorm(x=x, mean=meanlog, sd=sdlog), col="blue", lwd=2, add=T)
par(mfrow=c(1,1))

##################################################
# Compare the two methods of generating random deviates:
n <- 1000
niter <- 1000 # Resample 1000 times

x.mean.lnorm <- t(sapply(1:niter, function(i) {
		x1 <- rlnorm(n=n, meanlog=meanlog, sdlog=sdlog)
		x2 <- exp(meanlog + sdlog*rnorm(n=n, mean=0, sd=1))
		c(mean(x1), mean(x2))
		})
	)

x.sd.lnorm <- t(sapply(1:niter, function(i) {
		x1 <- rlnorm(n=n, meanlog=meanlog, sdlog=sdlog)
		x2 <- exp(meanlog + sdlog*rnorm(n=n, mean=0, sd=1))
		c(sd(x1), sd(x2))
		})
	)

apply(x.mean.lnorm,2,mean)
apply(x.mean.lnorm,2,sd)

apply(x.sd.lnorm,2,mean)
apply(x.sd.lnorm,2,sd)

dev.new(width=12, height=6)
par(mfrow=c(1,2))
	hist(x.mean.lnorm[,1], freq=FALSE, breaks=30, col="lightgrey")
	hist(x.mean.lnorm[,2], freq=FALSE, breaks=30, col="lightgrey")
par(mfrow=c(1,1))

##################################################
# Estimate moments of a lognormal distribution
# Requires
library(fitdistrplus)
##################################################

# Create random sample from lognormal distribution
n <- 10000
meanlog <- 1
sdlog <- 1/2

.x <- rlnorm(n=n, meanlog=meanlog, sdlog=sdlog)

plotdist(.x, breaks=30, col="red", type="l")

plotdist(log(.x), breaks=30, col="red", type="l")

descdist(.x, boot=999)
descdist(log(.x), boot=999)

f1.lnorm <- fitdist(.x, "lnorm") # Returns moments on log-scale
mean(.x); sd(.x) # For comparison

f1.norm <- fitdist(log(.x), "norm") # Same results as f1.lnorm!
mean(log(.x)); sd(log(.x)) # For comparison

summary(f1.lnorm)
summary(f1.norm)

gofstat(f1.lnorm)
gofstat(f1.lnorm)

plot(f1.lnorm)
plot(f1.norm)

# Bootstrap
b1 <- bootdist(f1.lnorm, niter=100)
summary(b1)

plot(b1)
points(median(b1$estim$meanlog), median(b1$estim$sdlog), pch=1, col="green", cex=2, lwd=4) # Bootstrap median
points(mean(log(.x)), sd(log(.x)), pch=16, col="red", cex=2) # Sample mean and sd
points(f1.lnorm$estimate[1],f1.lnorm$estimate[2],pch=15, col="blue", cex=2) # MLE
legend("topright", c("Bootstrap median","Sample estimates", "MLE estimates"), col=c("green", "red", "blue"), pch=c(1,16,15), bty="n")
##################################################
