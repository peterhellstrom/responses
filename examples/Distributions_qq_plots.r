###############################################################
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions/")
source("DistributionsFunctions.r")
###############################################################

n <- 1000
mu <- 4
sigma <- 3

x <- rnorm(n=n, mean=mu, sd=sigma)

dev.new()
par(mfrow=c(2,2))
	
	hist(x, breaks=30, col="steelblue", freq=FALSE)
	curve(dnorm(x, mu, sigma), add=TRUE, col=2, lwd=2)
	
	qqnorm(x, main="qqnorm")
	qqline(x, lty=2)
	title(sub="Line added with qqline")
	
	y <- rnorm(n,mean=mu,sd=sigma)
	
	qqplot(x,y,main="qqplot")
	abline(a=0, b=1, lty=2)
	title(sub="Line added with abline")
	
	qqplot(scale(x), y,main="qqplot, x standardized")
	qqline(x, lty=2)
	title(sub="Line added with qqline")
	
par(mfrow=c(1,1))

# Save output from qqplot
z <- qqplot(x, rnorm(n,mean=mu,sd=sigma), plot.it=FALSE)

###############################################################
# with GAMLSS (sample from a leptokurtic distribution)

n <- 1000
mu <- 4
sigma <- 3
nu <- 0.9

# Generate random variates
x <- rTF(n=n, mu=mu, sigma=sigma, nu=0.9)
# hist(x, col="steelblue", breaks=30)
range(x)

descdist(x, boot=999)

# Fit T-family distribution
f.TF <- gamlss(x~1, family="TF")

plot(f.TF) # Notice qq-plot in lower right corner

dev.new(width=12, height=6)
par(mfrow=c(1,2))
	qqnorm(scale(x)) # This is not what GAMLSS does
	qqline(scale(x), lty=2)
	
	qqnorm(resid(f.TF), col="darkgreen") # But instead, GAMLSS makes a qqplot of the residuals
	qqline(resid(f.TF), col="red")
par(mfrow=c(1,1))

# So if the GAMLSS-model provides a good fit, residuals should be normally distributed
dev.new()
hist(resid(f.TF), col="steelblue", breaks=30, freq=FALSE)
# Add standard normal pdf
curve(dnorm(x, mean=0, sd=1), col=2, lwd=2, add=TRUE, n=1001)

plot.new(); rqres.plot(f.TF, all=FALSE)
# Implementation in my custom function plot.Gamlss
dev.new(); plot.Gamlss(f.TF)

# Compare with normal distribution
f.NO <- gamlss(x~1, family="NO")

plot(f.NO)
plot.Gamlss(f.NO)

