##################################################
# Load library
library(fitdistrplus)
library(moments)
##################################################
# Generate a random sample from a theoretical distribution
# In this case, a normal distribution

n <- 50
x1 <- rnorm(n=n, mean=0, sd=1)

##################################################

# The package fitdistrplus uses qq-, pp-, density-, and cdf-plots extensively.
# It is therefore essential to understand these concepts!

# qq-plots are used for comparing a sample (measured distribution) to a reference distribution
# Check source code for qqnorm():
# https://svn.r-project.org/R/trunk/src/library/stats/R/qqnorm.R
# https://svn.r-project.org/R/trunk/src/library/stats/R/ppoints.R

# The data for a qq-plot is generated like this:

# Reference distribution = standard normal (qnorm)
qq.x <- qnorm(ppoints(n))[order(order(x1))]
qq.y <- x1

# ppoints generates a sequence of probability points:
# ppoints(n)
# qnorm gives the quantile function (the inverse of the cumulative density function)
# qnorm(ppoints(n))
# The double order-arguments makes sure that the quantile samples are matched against the sample.

# data.frame(qq.x, qq.y)

# Plot both data sets together, just to make sure we did it right
dev.new(width=12, height=6)
par(mfrow=c(1,2))
	plot(qq.x, qq.y)
	qqnorm(x1)
par(mfrow=c(1,1))

##################################################
# DISTRIBUTIONS IN R
# Remember how the different probablity distributions are defined in R:
# dnorm = density function (pdf)
# pnorm = distribution function, cumulative (cdf)
# qnorm = quantile function (inverse of the cdf)
# rnorm = random variates

# An example with the standard normal distribution
.x <- seq(-3,3,0.01)
dev.new()
par(mfrow=c(2,2))
	plot(.x, dnorm(.x), type="l", main="Density function"); abline(v=0, lty=2) # pdf
	plot(.x, pnorm(.x), type="l", main="Distribution function"); abline(v=0, lty=2, col=1) # cdf
	plot(ppoints(.x), qnorm(ppoints(.x)), type="l", main="Quantile function"); abline(h=0, lty=2) # qf
	hist(rnorm(length(.x)), breaks=30, col="lightgrey", main="Random variates", freq=FALSE); abline(v=0, lty=2) # random variates
	curve(dnorm(x), n=1001, add=T, lty=2)
par(mfrow=c(1,1))
##################################################

##################################################
# GRAPHICAL DISPLAY OF THE OBSERVED DISTRIBUTION

# Plot distribution

# The plotdist()-function by default plots a histogram and cumulative density function of the sample.
plotdist(x1)
# It is possible to change the usual plotting parameters:
plotdist(x1,col="red",type="b",pch=16)

# Plot a distrubution against data. 
# Note that the para argument must be a named list matching the arguments of the function chosen as distr.
# Using the function with distr and para defined return: histogram+density, qq-plot, cdf, pp-plot. 
plotdist(data=x1, distr="norm", para=list(mean=0, sd=1))

plotdist(data=x1, distr="norm", para=list(mean=0, sd=1), col="red", lwd=2, breaks=20) # Changing graphical parameters
# Plot against a different distribution:
plotdist(data=x1, distr="norm", para=list(mean=0, sd=2)) # Change the sd of the theoretical distribution

# The function plotdist produces different output for continuous and discrete data:
x2 <- rpois(n=50, lambda=1)
plotdist(x2, distr="pois", para=list(lambda=1))

plotdist(x2, discrete=FALSE)

##################################################
# CHARACTERIZATION OF THE OBSERVED DISTRIBUTION

# Descriptive statistics
# Here it is essential to remember skewness and kurtosis!
# Standard output is summary statistics and a Cullen-Frey graph.

descdist(x1, method="unbiased")
descdist(x1, method="sample") # Uses "common" functions for calculating skewness and kurtosis

# Kurtosis:
# Negative kurtosis - flat distribution - platykurtic
# Positive kurtosis - peaked distribution - leptokurtic
# Zero kurtosis - mesokurtic (i.e. normal)

# Skewness:
# Negative skewness - mean less than median - left skewed
# Positive skewness - mean larger than median - right skewed
# Zero skewness - mean equal to median - not skewed

kurtosis(x1)
skewness(x1)

# Kurtosis for a normal distribution is zero ONLY if it is corrected for excess kurtosis, i.e. subtracting 3:
kurtosis(rnorm(10000)) - 3

# Skewness and kurtosis are known to be non-robust.
# Apply non-parametric bootstrapping to take into account uncertainty in these estimators:
descdist(x1, boot=999)

##################################################
# FITTING OF A DISTRIBUTION

# These distributions are allowed:
# "norm", "lnorm", "pois", "exp", "gamma", "nbinom", "geom", "beta", "unif" and "logis"

# The standard method is to use maximum likelhood estimation, method="mle". Other options are available,
# see
?fdist

f1n <- fitdist(x1, "norm")
str(f1n)

summary(f1n) # Summary of fitted object, including parameter estimates, log-likelihood, and AIC.

# Goodness-of-fit, look for the smallest value. Can and will often reject the null hypothesis for large data sets, 
# mainly because of discrepancies in the tails of the distribution.

gofstat(f1n) # Computes goodness-of-fit statistics for a fit of a parametric distribution on non-censored data
gofstat(f1n, print=TRUE) # Prints p-values as well
plot(f1n) # Visually compare the sample with the theoretical distribution.

# Test using the moment matching estimator:
f1nb <- fitdist(x1, "norm", method="mme")
summary(f1nb)

# Try another distribution, uniform:
f1u <- fitdist(x1, "unif", method="mme")
summary(f1u)
gofstat(f1u)
plot(f1u)

##################################################
# SIMULATION OF THE UNCERTAINTY BY BOOTSTRAP
b1 <- bootdist(f1n)
summary(b1)
plot(b1)

b1p <- bootdist(f1n, method="param", niter=5001)
b1np <- bootdist(f1n, method="nonparm", niter=5001)

par(mfrow=c(1,2))
	plot(b1p)
	plot(b1np)
par(mfrow=c(1,1))

##################################################
##################################################
# Custom-defined or "non-basic" functions
library(NORMT3) # has the error function erf
library(fdrtool)

# Half-normal distribution
# pdf and cdf:
dhalfnorm1 <- function(x, theta) ((2*theta)/pi) * exp(-x^2 * theta^2 / pi)
phalfnorm1 <- function(q, theta) erf((theta*q) / sqrt(pi))
# For usage of my functions, it would be necessary to define the quantile function as qhalfnorm1

# The half-normal function is also available in package fdrtool

# Compare my definition of pdf and cdf with the functions in fdrtool:
dev.new(width=12, height=6)
par(mfrow=c(1,2))
curve(dhalfnorm(x=x, theta=2), from=0, to=2, n=1001)
curve(dhalfnorm1(x=x, theta=2), from=0, to=2, n=1001, lty=2, col=2, lwd=2, add=T)

curve(phalfnorm(q=x, theta=2), from=0, to=2, n=1001)
curve(phalfnorm1(q=x, theta=2), from=0, to=2, n=1001, lty=2, col=2, lwd=2, add=T)
par(mfrow=c(1,1))

# Simulate random variates from a half-normal distribution
n <- 150
.x <- rhalfnorm(n=n, theta=2)

# Plot a histogram of the sample
plotdist(.x, breaks=20, col="lightgrey", type="b")

# Compare the sample with different theoretical distributions;
# normal, exponential & half-normal
plotdist(.x, distr="norm", para=list(mean=mean(.x), sd=sd(.x)), breaks=30)

plotdist(.x, distr="exp", para=list(rate=1/mean(.x)), breaks=30)

plotdist(.x, distr="halfnorm", para=list(theta=1/mean(.x)), breaks=30)

# Fit different distributions
fx.halfnorm <- fitdist(.x, "halfnorm", start=list(theta=(1/mean(.x))))
fx.norm <- fitdist(.x, "norm")
fx.lnorm <- fitdist(.x, "lnorm")
fx.exp <- fitdist(.x, "exp")
fx.gamma <- fitdist(.x, "gamma")

# Create a list with fitted objects, and return AICc-table
fx.lst <- list(fx.halfnorm, fx.norm, fx.lnorm, fx.exp, fx.gamma)
extrDistr(fx.lst, order=TRUE)

# Plot some of the fitted objects
plot(fx.halfnorm, breaks=30)
plot(fx.gamma, breaks=30)
plot(fx.exp, breaks=30)