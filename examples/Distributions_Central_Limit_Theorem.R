library(fitdistrplus)

# Sample from a non-normal distribution, calculate mean for a large number of samples:

# Clearly non-normal variable, tri-modal
means <- c(10,15,20)
sds <- c(1,1,1)

.x <- rnorm(n=10000, mean=means, sd=sds)

# Do niter number of resamplings of size nsize from variable x, and calculate mean for each sample:
niter <- 10000
nsize <- 30

mean.x <- sapply(1:niter, function(i) mean(sample(x=.x, size=nsize, replace=TRUE)))

# Plot distribution of original variable x, and distribution of means of samples drawn from x
par(mfrow=c(1,2))
	hist(.x, freq=FALSE, breaks=30, main="Distribution of x", col="lightgrey")
	# for (i in 1:length(means)) curve(dnorm(x,means[i],sds[i]),add=T,n=1001,col=i, lwd=2)
	hist(mean.x, breaks=30, freq=FALSE, main=paste("Mean of", niter, "samples of size", nsize, "from x"), col="lightgrey")
	curve(dnorm(x,mean(mean.x),sd(mean.x)), n=1001, col=2, lty=1, lwd=2, add=T)
par(mfrow=c(1,1))

# Fit a normal distribution to the resamples:

f <- fitdist(mean.x, "norm")
dev.new()
plot(f, breaks=30)

# Sample from a normal distribution with mean=mu and sd=sigma
# Size of samples are 1:n, calculate mean for each sample.
# Purpose: demonstrate precision of sampling as sample size increases

# Define distribution parameters
mu <- 4
sigma <- 6

# Draw theoretical distribution
xlims <- c(floor(qnorm(0.025,mu,sigma)), ceiling(qnorm(0.975,mu,sigma)))

curve(dnorm(x,mu,sigma), xlim=xlims, ylab="Density")
lines(rep(mu,2), c(0,dnorm(mu,mu,sigma)), lty=2, col=2)

# Sample from the theoretical distribution
# Generate 1 sample for each sample size, 1:n
size <- 1:100
x <- lapply(size, function(i) rnorm(n=size[i], mean=mu, sd=sigma))

par(mfrow=c(1,2))
	# Mean
	plot(size, sapply(x, mean), xlab="n", ylab="mean")
	abline(h=mu, lty=2, col=2)
	# Add confidence limits
	lines(size, mu+(qnorm(0.025,0,1)*(sigma/sqrt(size))), lty=3, col=4)
	lines(size, mu+(qnorm(0.975,0,1)*(sigma/sqrt(size))), lty=3, col=4)
	# SD
	plot(size, sapply(x, sd), xlab="n", ylab="sd")
	abline(h=sigma, lty=2, col=2)
	lines(size, sigma+(qnorm(0.025,0,1)*(sigma/sqrt(size))), lty=3, col=4)
	lines(size, sigma+(qnorm(0.975,0,1)*(sigma/sqrt(size))), lty=3, col=4)
par(mfrow=c(1,1))

# Similar to above, but repeat nrepl number of time for each sample size i.

nrepl <- 100
size <- seq(10,10000,by=100)

x <- lapply(1:length(size), function(i) replicate(mean(rnorm(n=size[i], mean=mu, sd=sigma)), n=nrepl))

par(mfrow=c(1,2))
	plot(size, sapply(x, mean), xlab="n", ylab="mean")
	abline(h=mu, lty=2, col=2)

	plot(size, sapply(x, sd), xlab="n", ylab="sd")
	lines(size, sigma/sqrt(size), lty=3, col=4)
par(mfrow=c(1,1))


# Demo of Central Limit Theorem 3
# Show the connection with sample size and sample accuracy.

# Sample from normal distribution
n <- c(1,4,16,64,100)
mu <- 5
sigma <- 16
nrepl <- 1000


d <- t(replicate(sapply(1:length(n), function(i) mean(rnorm(n[i], mean=mu, sd=sigma))), n=nrepl))
colnames(d) <- n

apply(d,2,mean)
apply(d,2,sd)
sigma/sqrt(n)

d.ls <- lapply(1:length(n), function(i) density(d[,i]))
names(d.ls) <- n

ylims <- range(sapply(1:length(n), function(i) range(d.ls[[i]]$y)))

plot(x=range(d), y=ylims, type="n", xlab="x", ylab="Density", font.lab=1, las=1, main="Central limit theorem")

for (i in 1:length(n)) {
	lines(d.ls[[i]],col=i ,lty=2)
	curve(dnorm(x, mean=mu, sd=sigma/sqrt(n[i])),add=T,col=i)
}


# Sample from exponential distribution

n <- c(1,4,16,64,100)
rate <- 1/4
nrepl <- 100

d <- t(replicate(sapply(1:length(n), function(i) mean(rexp(n[i], rate=rate))), n=nrepl))
colnames(d) <- n

apply(d,2,mean)
apply(d,2,sd)
(1/rate)/sqrt(n)

d.ls <- lapply(1:length(n), function(i) density(d[,i]))
names(d.ls) <- n

ylims <- range(sapply(1:length(n), function(i) range(d.ls[[i]]$y)))

plot(x=range(d), y=ylims, type="n", xlab="x", ylab="Density", font.lab=1, las=1, main="Central limit theorem")

for (i in 1:length(n)) {
	lines(d.ls[[i]],col=i ,lty=2)
	curve(dnorm(x, mean=1/rate, sd=(1/rate)/sqrt(n[i])), add=TRUE, col=i)
}
