###############################################################
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions")
source("DistributionsFunctions.r")
###############################################################

# Define parameters for 3-parameter t-distribution
mu <- 5
sigma <- 2
nu <- 1.5
nu <- 1 # if nu=1, then t = cauchy

# If you set nu=1, you'll se that the cauchy distribution is a special case of the t-family distribution!
###############################################################
# Probability density function
#qs <- c(0.025, 0.975)
qs <- c(0.05, 0.95)
xlims <- qcauchy(qs,mu,sigma)

curve(dTF(x, mu=mu, sigma=sigma, nu=nu), xlim=xlims, n=1001, lwd=2, ylab="Density", main="3-parameter t- & Cauchy distribution")
curve(dcauchy(x, location=mu, scale=sigma), n=1001, col=2, lwd=2, lty=2, add=T)
lines(x=c(mu,mu), y=c(0,dTF(x=mu, mu=mu,sigma=sigma,nu=nu)),col=4,lty=3)
legend("topright", c("t","Cauchy"), col=c(1,2), lwd=c(2,2), lty=c(1,2), bty="n", title="Distribution")
# Add mathematical notation to plot
#text(-5,0.15, bquote(paste(mu== .(mu), ", ", sigma==.(sigma), ", ", nu==.(nu))))
legend("topleft", legend=bquote(paste(mu== .(mu), ", ", sigma==.(sigma), ", ", nu==.(nu))), bty="n")

###############################################################
# GAMLSS does not have a specified Cauchy-family, but it should nevertheless be possible to fit
# a cauchy-distribution to data by constraining the parameter nu.

# Fit untruncated distributions

# Generate data
n <- 2500
x <- rcauchy(n=n, location=mu, scale=sigma)

# Fit Cauchy-distribution with fitdistrplus
fm1 <- fitdist(x, "cauchy")
summary(fm1)
#plot(fm1, breaks=30, col="steelblue")

# Fit Cauchy-distribution with gamlss, constrain parameter nu to 1.
fm2 <- gamlss(x~1, family="TF", nu.start=1, nu.fix=TRUE)
summary(fm2)

# Fit TF-distribution with gamlss, unconstrained
fm3 <- gamlss(x~1, family="TF")
summary(fm3)

# Inspect parameter estimates:
fm1$estimate
c(location=fm2$mu.coefficients, scale=exp(fm2$sigma.coefficients))
extrcoef.gamlss(fm3)

###############################################################

# Fit truncated distributions

# Generate truncated data
n <- 2500
a <- 0
b <- 1300

x <- rtrunc(n=n, spec="cauchy", a=a, b=b, location=mu, scale=sigma)
range(x)

# Fit Cauchy-distribution with fitdistrplus
fm1.tr <- fitdist(x, "cauchy")
summary(fm1.tr)
#plot(fm1.tr, breaks=30, col="steelblue")

# Fit left- & right-truncated Cauchy-distribution with gamlss, constrain parameter nu to 1.
gen.trun(par=c(a,b), family="TF", type="both")
fm2.tr <- gamlss(x~1, family=trun(par=c(a,b), family="TF", type="both"), nu.start=1, nu.fix=TRUE)
summary(fm2.tr)

# Fit left- & right-truncated TF-distribution with gamlss, unconstrained
gen.trun(par=c(a,b), family="TF", type="both")
fm3.tr <- gamlss(x~1, family=trun(par=c(a,b), family="TF", type="both"))
summary(fm3.tr)

# Fit right-truncated TF distribution
gen.trun(par=b, family="TF", type="right")
fm4.tr <- gamlss(x~1, family=trun(par=b, family="TF", type="right"))
summary(fm4.tr)

# Inspect parameter estimates:
fm1.tr$estimate
c(location=fm2.tr$mu.coefficients, scale=exp(fm2.tr$sigma.coefficients))
extrcoef.gamlss(fm3.tr)
extrcoef.gamlss(fm4.tr)

hist(x, breaks=30, col="steelblue",freq=FALSE)
lines(density(x))
curve(dtrunc(x, spec="cauchy", a=a, b=b, location=mu, scale=sigma),add=T,n=1001,col=2,lwd=2)

###############################################################

