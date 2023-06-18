# Random variates generation
# rnorm vs. sample

n <- 10000
nrepl <- 1000

mu <- 0
sigma <- 1

system.time(x1 <- replicate(rnorm(n=n,mu,sigma),n=nrepl))

fac <- 100
system.time({
	xs <- rnorm(n=n*fac, mu, sigma)
	x2 <- replicate(sample(xs,size=n,replace=TRUE), n=nrepl) })

x1.m <- colMeans(x1)	
x2.m <- colMeans(x2)

# bias

plot(density(x1.m - mu))
lines(density(x2.m - mu), col=2)
abline(v=mean(x1.m - mu), col=1, lty=2)
abline(v=mean(x2.m - mu), col=2, lty=2)

mean(x1.m - mu)
mean(x2.m - mu)

####
# The additive law of expectation
# http://en.wikipedia.org/wiki/Sum_of_normally_distributed_random_variables
# Error accumulation
# Distribution of sums of random variates
# Take nperobj random variates per obj from a normal distribution
# Replicate this for

mu <- 20
sigma <- 5

nobj <- 20000
nperobj <- 6

# Approach 1, generate nobj * nperobj random variates
erracc <- function(nobj, nperobj, mean, sd) replicate(sum(rnorm(n=nperobj, mean=mean, sd=sd)), n=nobj)

system.time(test1 <- erracc(nobj=nobj, nperobj=nperobj, mean=mu, sd=sigma))

# Approach 2 - shortcut. Generate nobj random variates
# So instead of taking nperobj observation, isn't it possible to sample from the expected
# distribution of sums of normally distributed values?
# 40-50 times faster.

mu.sum <- sum(rep(mu,nperobj)) # expected mean of sum
sigma.sum <- sqrt(sum(rep(sigma^2,nperobj))) # expected sd of sum

system.time(test2 <- rnorm(n=nobj, mu.sum, sigma.sum))

xlims <- c(qnorm(0.0001, mu.sum, sigma.sum), qnorm(0.9999, mu.sum, sigma.sum))
curve(dnorm(x, mu.sum, sigma.sum), col=2, xlim=xlims, lwd=2, main=expression(paste("Same ", mu, " & ", sigma))) # Expectation
lines(density(test1),col=3) # Simulation, method 1
lines(density(test2), col=4) # Simulation, method 2 (shortcut)
legend("topleft", c("Expectation","Simulation, Long","Simulation, Shortcut"), lty=c(1,1,1), lwd=c(2,1,1), col=c(2,3,4), bty="n", cex=0.8)


# if means and sds differ, the same theorem applies:
nperobj <- 3
mus <- c(0,-1,3)
sigmas <- c(1,2,1)

if(length(mus) != length(sigmas) | length(mus) != nperobj) stop("TEST")

mu.sums <- sum(mus)
sigma.sums <- sqrt(sum(sigmas^2))

test3 <- erracc(nobj=nobj, nperobj=nperobj, mean=mus, sd=sigmas)
test4 <- rnorm(n=nobj, mu.sums, sigma.sums)

xlims <- c(qnorm(0.0001, mu.sums, sigma.sums), qnorm(0.9999, mu.sums, sigma.sums))
curve(dnorm(x, mu.sums, sigma.sums), col=2, xlim=xlims, lwd=2, main=expression(paste("Different ", mu, " & ", sigma))) # Expectation
lines(density(test3),col=3) # Simulation, method 1
lines(density(test4), col=4) # Simulation, method 2 (shortcut)
legend("topleft", c("Expectation","Simulation, Long","Simulation, Shortcut"), lty=c(1,1,1), lwd=c(2,1,1), col=c(2,3,4), bty="n", cex=0.8)

