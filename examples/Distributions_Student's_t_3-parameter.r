###############################################################
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions/")
source("DistributionsFunctions.r")
###############################################################

# Generate random variates from the three-parameter t-distribution.
# family = "TF" in gamlss

n <- 1000

mu <- 4.667
sigma <- exp(-0.8761)
nu <- exp(0.4338)

xtf <- rTF(n=n, mu=mu, sigma=sigma, nu=nu)
range(xtf)
hist(xtf, breaks=200, col="steelblue", freq=FALSE, xlim=range(xtf)/2)
curve(dTF(x, mu=mu, sigma=sigma, nu=nu),add=T,col=2,lwd=2)

# Truncated distribution
# Truncation points
a <- 0.69
b <- 9.5
trpars <- c(a,b)

xtf.tr <- rtrunc(n=n, spec="TF", mu=mu, sigma=sigma, nu=nu, a=a, b=b)
plotdist(xtf.tr, breaks=50, col="steelblue")

range(exp(xtf.tr))
quantile(exp(xtf.tr))
IQR(exp(xtf.tr))
hist(exp(xtf.tr), breaks=50, col="steelblue", freq=FALSE)

###############################################################
# Probability density function (custom function dt3, see DistributionsFunctions.r)
xlims <- c(mu-5, mu+5)
curve(dTF(x, mu=mu, sigma=sigma, nu=nu), xlim=xlims, n=1001, lwd=2, ylab="Density", main="3-parameter t-distribution")
curve(dt3(x, mu=mu, sigma=sigma, nu=nu), col=2, lwd=2, lty=4, add=T, n=1001)
lines(density(xtf), col=3)

###############################################################
# HOW TO generate random variates from the 3-parameter t-distribution in gamlss with R's standard rt() function
# y ~ t(mu,sigma,nu) == y ~ tstand(nu)*sigma + mu
# Random variates (custom function rt3tr, see DistributionsFunctions.r)

n <- 10000
mu <- 4.667
sigma <- exp(-0.8761)
nu <- exp(0.4338)
a <- 0.69; b <- 9.5
trpars <- c(a,b)

xtf.tr <- rt3tr(n=n, a=a, b=b, mu=mu, sigma=sigma, nu=nu)
range(xtf.tr) # Should ~equal trpars

# Plot histogram of random variates
hist(xtf.tr, breaks=100, col="steelblue", freq=FALSE, xlim=c(-1,10))
# Plot untruncated t-family distribution
curve(dTF(x, mu=mu, sigma=sigma, nu=nu),add=T, col=2, lwd=2)
# Plot truncated t-family distribution
gen.trun(trpars, family="TF", type="both") # generate distribution functions
curve(dTFtr(x, mu=mu, sigma=sigma, nu=nu), xlim=trpars, add=T, col=3, lwd=2)

# Random variates with GAMLSS
xtf.tr.gamlss <- rTFtr(n=n, mu=mu, sigma=sigma, nu=nu)
plotdist(xtf.tr.gamlss, breaks=50, col="steelblue")

###############################################################
# EXTRAS
###############################################################
# CHECK consistency betwwen rt3tr & rTFtr
n <- 1000
nrepl <- 10000

xtf.test <- replicate({
	xtf.tr <- rt3tr(n=n, a=a, b=b, mu=mu, sigma=sigma, nu=nu)
	xtf.tr.gamlss <- rTFtr(n=n, mu=mu, sigma=sigma, nu=nu)
	mean(xtf.tr) - mean(xtf.tr.gamlss)
	}, n=nrepl)
	
hist(xtf.test, breaks=50, col="steelblue", freq=FALSE)
lines(density(xtf.test), col=2, lwd=2)

###############################################################
# Check computation time for the different functions

n <- 1000
nrepl <- 100

system.time(d0a <- rt(n=n, df=nu))
system.time(d0b <- (rt(n=n, df=nu) * sigma) + mu)
system.time(d0c <- rTF(n=n, mu=mu, sigma=sigma, nu=nu))

system.time(d1a <- (rtrunc(n=n, spec="t", df=nu, a=a, b=b) * sigma) + mu)
system.time(d1b <- rt3tr(n=n, a=a, b=b, mu=mu, sigma=sigma, nu=nu))
system.time(d1c <- rTFtr(n=n, mu=mu, sigma=sigma, nu=nu))

system.time(d1r <- replicate(rt3tr(n=n, a=a, b=b, mu=mu, sigma=sigma, nu=nu), n=nrepl))
system.time(d2r <- replicate(rTFtr(n=n, mu=mu, sigma=sigma, nu=nu), n=nrepl))

###############################################################
