# Monte Carlo methods are based on generation of random variates.
# But how should the random variates be generated?

# Generate nrepl number of unique data sets and calculate mean & sd for each unique set (method="random")
# OR
# Generate one (big) dataset and sample with replacement from that data set (method="sample")

# Both the law of large numbers and the central limit theorem applies:
# If n is large, the true mean should be accurately estimated
# However, if n is small, the resampling in method="sample" will estimate the mean of the single sample.
# According to the clt, sd of a random sample from a normal distribution is ~ sigma/sqrt(n)
# The method="sample" is therefore biased for 'small' samples. One way to remedy this shortcoming is
# to generate a very large sample with random variates (n = nlarge), and then sample with replacement from the larger sample (size = n << nlarge)
# These three appraoches are outlined below:

clt.norm <- function(n=100, nlarge=NULL, nrepl=10, mu=0, sigma=1, method="random") {
	
	if (method == "random") {
		xr <- t(replicate( {
			x <- rnorm(n=n, mean=mu, sd=sigma)
			c(mean=mean(x), sd=sd(x))
		}, n=nrepl))
	}

	if (method == "sample") {
		if (is.null(nlarge)) nlarge=n
		x <- rnorm(n=nlarge, mean=mu, sd=sigma)
		xr <- t(replicate({
			x.s <- sample(x, size=n, replace=TRUE)
			c(mean=mean(x.s), sd=sd(x.s))
		}, n=nrepl))
	}
	data.frame(xr)
}

n <- 1000
nrepl <- 10000
nlarge <- 20000
mu <- 0
sigma <- 1


system.time(test1 <- clt.norm(n=n, nrepl=nrepl, mu=mu, sigma=sigma, method="random"))
system.time(test2 <- clt.norm(n=n, nrepl=nrepl, mu=mu, sigma=sigma, method="sample"))
system.time(test3 <- clt.norm(n=n, nrepl=nrepl, nlarge=nlarge, mu=mu, sigma=sigma, method="sample"))

d1 <- density(test1$mean)
d2 <- density(test2$mean)
d3 <- density(test3$mean)

curve(dnorm(x, mean=mu, sd=sigma/sqrt(n)), n=1001, xlim=c(-.15,.15), lty=2, lwd=2, col="steelblue")
lines(d1, col=1)
lines(d2, col=2)
lines(d3, col=3)
legend("topleft", c("Theoretical", "Random", "Sample, n", "Sample, nlarge"), col=c("steelblue",1,2,3), lwd=c(2,1,1,1), bty="n")
# Difference between d2 and d3 is controlled by the ratio nlarge/n



graphics.off()