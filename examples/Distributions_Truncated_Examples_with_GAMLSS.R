###############################################################
source("C:/WORK/ANALYSER/-= R Development =-/distRibutions/distRibutions.r")
###############################################################
# Examples from this paper (where the code was published):
# Nadarajah & Kotz 2006 j stat soft

x <- seq(-3, 3, by = 0.1)
## standard normal pdf
y1 <- dnorm(x)
## truncated standard normal pdf at a=-0.5 & b=0.5
y2 <- dtrunc(x, "norm", a = -0.5, b = 0.5, mean = 0, sd = 2)
## truncated standard normal pdf at a=-1 & b=1
y3 <- dtrunc(x, "norm", a = -1, b = 1, mean = 0, sd = 2)
## truncated standard normal pdf at a=-2 & b=2
y4 <- dtrunc(x, "norm", a = -2, b = 2, mean = 0, sd = 2)
yrange <- range(y1, y2, y3, y4)
plot(x, y1, type = "l", xlab = "x", ylab = "PDF", xlim = c(-3,3), ylim = yrange)
lines(x, y2, lty = 2)
lines(x, y3, lty = 3)
lines(x, y4, lty = 4)

## generate 100 random numbers from the truncated Weibull with a=1, b=2
set.seed(1)
x <- rtrunc(100, "weibull", a = 1, b = 2, shape = 2)
## compute the observed order statistics
x <- sort(x)
## compute the corresponding expected values
y <- qtrunc((seq(1,100) - 0.375)/(100 + 0.25), "weibull", a = 1, b = 2, shape = 2)
m1 <- min(x, y)
m2 <- max(x, y)
plot(x, y, xlim = c(m1, m2), ylim = c(m1, m2), xlab = "Observed", ylab = "Expected")
abline(0, 1)

###############################################################
###############################################################
# My own examples:
###############################################################

####
# Normal distribution

mu <- 0; sigma <- 1
xlims <- c(-3,3); ylims <- c(0,0.65)

# Truncation points
a <- c(-3,-2,-1); b <- c(3,2,1)

# Test assymmetric truncation points
a <- c(-3,-2,-1); b <- c(2,1.5,0.75)

curve(dnorm(x, mean=mu, sd=sigma), xlim=xlims, ylim=ylims, n=1001, ylab="Density", main="Truncated normal distribution")
for (i in 1:length(a)) curve(dtrunc(x, "norm", a=a[i], b=b[i], mean=mu, sd=sigma), xlim=xlims, add=T, lty=(i+1), col=(i+1), n=1001)


# Compare halfnormal with normal left-truncated at zero
mu <- 0; sigma <- 1
xlims <- c(0,3); ylims <- c(0,0.85)
a <- c(0,0,0); b <- c(2,3,4)

curve(dhalfnorm(x, theta=sqrt(pi/2)), xlim=xlims, ylim=ylims, n=1001, ylab="Density", main="Truncated normal distribution")
for (i in 1:length(a)) curve(dtrunc(x, "norm", a=a[i], b=b[i], mean=mu, sd=sigma), xlim=xlims, add=T, lty=(i+1), col=(i+1), n=1001)

####
# Lognormal distribution
# a and b NOT on log-scale, but meanlog and sdlog are!

mu <- 2
sigma <- 0.6
xlims <- c(0,100)
ylims <- c(0,0.2)

# Truncation points
a <- c(0,5,25)
b <- c(15,30,500)

curve(dlnorm(x, meanlog=mu, sdlog=sigma), xlim=xlims, ylim=ylims, n=1001, ylab="Density", main="Truncated lognormal distribution")
for (i in 1:length(a)) curve(dtrunc(x, "lnorm", a=a[i], b=b[i], meanlog=mu, sdlog=sigma), xlim=xlims, add=T, lty=(i+1), col=(i+1), n=1001)

xr <- rtrunc(n=1000, spec="lnorm", a=25, b=500, meanlog=2, sdlog=0.6)
mean(xr)
quantile(xr, c(0.025, 0.5, 0.975))
plot(fitdist(xr, "lnorm"), breaks=30, col="steelblue") # untruncated distribution

####
# Generate random variates, and compare with theoretical quantile distribution
n <- 1000
a <- 0.25
b <- 20
meanlog <- 0.5
sdlog <- 1.08

# Generate random variates from truncated distribution
.x <- rtrunc(n=n, "lnorm", a=a, b=b, meanlog=meanlog, sdlog=sdlog)
range(.x)
mean(log(.x)); sd(log(.x))

# Quantile distribution
.x.E <- qtrunc(p=ppoints(n), spec="lnorm", a=a, b=b, meanlog=meanlog, sdlog=sdlog) # Truncated lognormal
.x.E.noT <- qlnorm(p=ppoints(n), meanlog=meanlog, sdlog=sdlog) # Untruncated lognormal

dev.new()
par(mfrow=c(2,2))
	# Plot histogram of random variates
	hist(.x, breaks=30, freq=FALSE, col="steelblue", main="x")
	lines(density(.x))
	# Plot histogram of log(random variates)
	hist(log(.x), breaks=30, freq=FALSE, col="steelblue", main="log(x)")
	lines(density(log(.x)))
	abline(v=meanlog, lwd=2, lty=2, col=2)
	# Plot ecdf
	plot(ecdf(.x), main="Empirical cdf")
	# Plot qq-plot
	plot(sort(.x),.x.E, type="p", pch=".", cex=3, main="qq-plot") # truncated q
	points(sort(.x), .x.E.noT, col=2, type="p", pch="+", cex=0.5) # untruncated q
	abline(0,1,lty=2,lwd=1)
	legend("bottomright", c("Truncated","Untruncated"), pch=c(".","+"), col=c(1,2), cex=1.3, bty="n")
par(mfrow=c(1,1))

###############################################################
# Estimate shape and scale parameters for a truncated distribution.
# Use the gamlss.tr package

trpoints <- c(0.25, 20)
n <- 10000
mu <- 0.5; sigma <- 1.08
# Generate functions for a truncated distribution
gen.trun(par=trpoints, family="LOGNO", type="both")
# Create a sample from the generated distribution
sam <- rLOGNOtr(n, mu=mu, sigma=sigma)

mNO <- histDist(log(sam), family="NO", density=TRUE)
plot(mNO)
wp(mNO)
# Q.stats(mNO)

mean(sam); sd(sam)
mean(log(sam)); sd(log(sam))

# Not possible to use truncated distributions with fitdistrplus
# Not compatible with gamlss.tr:
# fitdist(sam, distr="LOGNOtr", start=list(mu=mu,sigma=sigma)) # does NOT work, since the created function LOGNOtr has no arguments, see args(LOGNOtr)

dev.new(width=12, height=6)
par(mfrow=c(1,2))
	hist(sam, breaks=30, col="steelblue", freq=FALSE)
	#curve(dtrunc(x, spec="lnorm", a=trpoints[1], b=trpoints[2], meanlog=mu, sdlog=sigma), add=TRUE, lwd=2, col=2, n=1001)
	curve(dLOGNOtr(x,mu=mu,sigma=sigma), xlim=trpoints, n=1001, lwd=2, col=2, add=T)
	hist(log(sam), breaks=30, col="steelblue", freq=FALSE)
	curve(dtrunc(x, spec="norm", a=log(trpoints[1]), b=log(trpoints[2]), mean=mu, sd=sigma), add=TRUE, lwd=2, col=2, n=1001)
par(mfrow=c(1,1))

# fit the distribution to the data
mod1 <- gamlss(formula=sam~1, sigma.formula=~1, family=trun(par=trpoints, family="LOGNO", type="both"))
mod1
summary(mod1); coef(mod1)
plot(mod1)

# Extract fitted parameters (note that link function for sigma is log!)
fpars <- as.numeric(c(mod1$mu.coefficients, exp(mod1$sigma.coefficients)))
fpars

# Effect of choosing a higher right-truncation point?
trpoints1a <- c(0.25, 25)
gen.trun(par=trpoints1a, family="LOGNO", type="both") # Must be defined again!!!
mod1a <- gamlss(formula=sam~1, sigma.formula=~1, family=trun(par=trpoints1a, family="LOGNO", type="both"))
summary(mod1a)
coef(mod1); coef(mod1a)
plot(mod1a)

fpars1a <- as.numeric(c(mod1a$mu.coefficients, exp(mod1a$sigma.coefficients)))
fpars; fpars1a


plot(density(sam), ylim=c(0,0.45))
curve(dLOGNOtr(x,mu=fpars[1],sigma=fpars[2]),xlim=trpoints,n=1001,add=T,lty=2,col=2)
curve(dLOGNOtr(x,mu=fpars1a[1],sigma=fpars1a[2]),xlim=trpoints1a,n=1001,add=T,lty=2,col=4)
legend("topright", c("density", paste("right-truncated at", trpoints[2]), paste("right-truncated at", trpoints1a[2])),
	col=c(1,2,4), lty=c(1,2,2), bty="n")

# now create a gamlss.family object before the fitting for comparison
LOGNOtrunc.Both <- trun(par=trpoints, family="LOGNO", type="both", local=FALSE)
mod2 <- gamlss(sam~1, family=LOGNOtrunc.Both)

summary(mod2)
coef(mod2)

###############################################################

# May be easier to use functions in the file truncated.r, rather than using gamlss.tr...
# This code validates/compares the d and r functions

# 1) Compare prob density function, d
# Matches perfectly:
mu <- 3
sigma <- 1.08
fpars <- c(mu, sigma)
trpars <- c(a=0.25, b=20.75)

gen.trun(par=trpars, family="LOGNO", type="both")

dev.new()
curve(dLOGNOtr(x, mu=fpars[1], sigma=fpars[2]), xlim=trpars, n=1001, lwd=2)
curve(dtrunc(x, spec="lnorm", a=trpars["a"], b=trpars["b"], meanlog=fpars[1], sdlog=fpars[2]), 
	from=trpars["a"], to=trpars["b"], n=1001, col=2, lty=2, add=T)

# 2) Compare random variates generation, r
n <- 1000
nrepl <- 10000

# Calculate difference between data sets generated´by rtrunc and rLOGNOtr
xdiff <- sapply(1:nrepl, function(i) {
	x1 <- rtrunc(n, spec="lnorm",a=trpars["a"], b=trpars["b"], meanlog=fpars[1], sdlog=fpars[2])
	x2 <- rLOGNOtr(n=n, mu=fpars[1], sigma=fpars[2])
	diffx1x2 <- mean(x1) - mean(x2)
	diffx1x2
	})

dev.new() # should be symmetrically centered around 0
hist(xdiff, freq=FALSE, breaks=30, col="steelblue", xlab="diff(x1-x2)", main=paste("mean difference (x1-x2) = ", round(mean(xdiff),4), sep=""))
title(sub=paste(nrepl, "simulations, each with sample size n =", n))
lines(density(xdiff))
abline(v=mean(xdiff), lty=2, col=2, lwd=2)

###############################################################
