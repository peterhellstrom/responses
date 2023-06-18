# Adapted from: http://www.r-bloggers.com/bootstrapping-the-truncated-normal-distribution/

#### Truncated Normal distribution ####
# & Inverse Mills ratio
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions")
source("DistributionsFunctions.r")

mu <- 0
sigma <- 4
tr.left <- -2 # left truncated at this x-value

# Plot pdf
xlims <- c(-16,16)
curve(dtrunc(x, spec="norm", a=tr.left, mean=mu, sd=sigma), xlim=xlims, col=2, n=1001, 
	main="Truncated Normal Distribution", xlab="x", ylab="Density")
curve(dnorm(x, mean=mu, sd=sigma), add=T, col=1, n=1001)
abline(v=tr.left, col=1, lty=2 ,lwd=2)
legend("topleft",legend=c("Original","Truncated","Truncation\nPoint"),col=c(1,2,1),lty=c(1,1,2),lwd=c(1,1,2), bty="n")

# Plot cdf
curve(ptrunc(x, spec="norm", a=tr.left, mean=mu, sd=sigma), xlim=xlims, col=2, n=1001,
	main="Truncated Normal Distribution", xlab="x", ylab="Cumulative density")
curve(pnorm(x, mean=mu, sd=sigma), add=T, col=1, n=1001)
abline(v=tr.left, col=1, lty=2 ,lwd=2)
legend("topleft",legend=c("Original","Truncated","Truncation\nPoint"),col=c(1,2,1),lty=c(1,1,2),lwd=c(1,1,2), bty="n")

# Create 'nbig" number of random normal variates and take repeated samples of size n, replicate nrepl number of times.
nbig <- 100000
nrepl <- 10000
n <- 200

x.source <- rnorm(n=nbig, mean=mu, sd=sigma)
s.sample <- replicate(sample(x.source, size=n, replace=TRUE), n=nrepl)

x.source.tr <- rtrunc(n=nbig, spec="norm", a=tr.left, mean=mu, sd=sigma)
s.sample.tr <- replicate(sample(x.source.tr, size=n, replace=TRUE), n=nrepl)

# trunc.sum is our row means for various outputs from the bootstrap
trunc.sum <- matrix(0,nrepl,3)
trunc.sum[,1] <- colMeans(s.sample)
trunc.sum[,2] <- colMeans(s.sample.tr)

# Alpha is a scale factor. It gives us a scale free measure of the truncation point
alpha <- (tr.left - mu) / sigma
# Calculate the Inverse Mills ratio
lambda1 <- dnorm(alpha,mean=0,sd=1) / (1 - pnorm(alpha, mean=0,sd=1))
lambda2 <- -dnorm(alpha,mean=0,sd=1) / (pnorm(alpha, mean=0,sd=1))

mu + sigma*lambda1
mu + sigma*lambda2

####
xlims <- range(range(density(s.sample)$x), range(density(s.sample.tr)$x))
ylims <- range(range(density(s.sample)$y), range(density(s.sample.tr)$y))

plot(x=xlims, y=c(0,ylims[2]+0.05), type="n",main=paste(nrepl, "Random Samples"), xlab="x", ylab="Density")
lines(density(s.sample.tr), col=2)
lines(density(s.sample), col=1)
curve(dtrunc(x, spec="norm", a=tr.left, mean=mu, sd=sigma), n=1001, add=T, col=2, lty=2)
curve(dnorm(x, mean=mu, sd=sigma), n=1001, add=T, col=1, lty=2)
abline(v=mean(s.sample.tr), col=4, lty=2)
abline(v=mu + sigma*lambda1, col=3, lty=3)

#for (i in 1:nrepl) lines(density(s.sample.tr[,i]), lwd=0.2, col=2)
#for (i in 1:nrepl) lines(density(s.sample[,i]), lwd=0.2, col=1)

####
hist(trunc.sum[,2],breaks=50, col="steelblue", freq=FALSE, main="Mean of E(x)")
abline(v=mu + sigma*lambda1, lwd=2, col=2) # Inverse Mills Ratio, expected value for x
legend("topleft", c("Inverse Mills Ratio"), lwd=2, col=2, bty="n")

mean(trunc.sum[,2])
mu + sigma*lambda1

####
