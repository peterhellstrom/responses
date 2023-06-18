##################################################
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions")
source("DistributionsFunctions.r")
##################################################

# Pareto distribution type IV
# Burr distribution a special case of the Pareto type IV distribution.
# EasyFit parameters: k=shape, alpha=1/inequality, beta=scale, gamma=location

# Example
loc <- 0
sc <- 4.6189
ineq <- 1/12.248
shp <- 0.99908

# The Burr distribution is available in:
# VGAM {d,p,q,r}paretoIV, with location=0
# gamlss, as a special case of the Generalized Beta type 2, with nu set to 1.

# Draw a Burr distribution, and add a truncated version:
curve(dparetoIV(x, location=loc, scale=sc, inequality=ineq, shape=shp), from=0, to=9, n=1001, ylab="Density", ylim=c(0,0.7))
curve(dtrunc(x, spec="paretoIV", a=3, b=6, location=loc, scale=sc, inequality=ineq, shape=shp), from=0, to=9, n=1001, add=T, col=2, lty=2)

# Check that the ParetoIV and Generalized Beta 2 parameterizations are indeed equivivalent:
curve(dparetoIV(x, location=loc, scale=sc, inequality=ineq, shape=shp), from=0, to=9, n=1001, ylab="Density", ylim=c(0,0.7))
curve(dGB2(x, mu=sc, sigma=1/ineq, nu=1, tau=shp), add=T, col=2, n=1001, lty=2)

# Generate data from a Burr distribution
n <- 10000
x <- rparetoIV(n=n, location=loc, scale=sc, inequality=ineq, shape=shp)
range(x)

hist(x, breaks=30, col="steelblue", freq=FALSE)
lines(density(x))

# Test 
# x2 <- rtrunc(n, spec="paretoIV", a=0, b=100, location=loc, scale=sc, inequality=ineq, shape=shp)
# lines(density(x2), col=2)

# Fit a model with gamlss
fm1 <- gamlss(x~1, family="GB2", nu.fix=TRUE, nu.start=1, n.cyc=500)
# summary(fm1)

# ?GB2
exp(fm1$mu.coefficients) # scale
fm1$sigma.coefficients # 1/inequality
exp(fm1$nu.coefficients) # 1
exp(fm1$tau.coefficients) # shape

##################################################
# Fit a truncated Burr-distribution
# Not quite sure about this, mu estimates are totally off...

# Generate data
n <- 1000
a <- 2
b <- 6
trpars <- c(a,b)

xt <- rtrunc(n, spec="paretoIV", a=a, b=b, location=loc, scale=sc, inequality=ineq, shape=shp)

plot(density(xt))
range(xt)

gen.trun(par=trpars, family="GB2", type="both")

gen.trun(par=trpars, family="GB2", type="both", nu.fix=TRUE, nu.start=1)

fm2 <- gamlss(xt~1, family=trun(par=trpars, family="GB2", type="both"), nu.fix=TRUE, nu.start=1)
summary(fm2)

exp(fm2$mu.coefficients) # scale
fm2$sigma.coefficients # 1/inequality
exp(fm2$nu.coefficients) # 1
exp(fm2$tau.coefficients) # shape

hist(xt, breaks=30, col="steelblue", freq=FALSE)
lines(density(xt))
curve(dtrunc(x, spec="paretoIV", a=a, b=b, location=loc, scale=sc, inequality=ineq, shape=shp), add=T, col=2)

curve(dGB2tr(x, mu=exp(fm2$mu.coefficients), sigma=fm2$sigma.coefficients, nu=1, tau=exp(fm2$tau.coefficients)), add=T, col=4)
