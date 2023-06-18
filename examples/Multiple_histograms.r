# Plot multiple histograms

library(plotrix)
library(MASS)

# Generate data
n <- 200 # sample size
x <- rnorm(n=n, 0, 1) # random variates
g <- factor(rep(1:2, each=100)) # grouping variable

xnew <- split(x,g) # split values in x, conditional on grouping variable

# plotrix
dev.new()
par(mfrow=c(2,2))
	mh <- multhist(xnew, beside=FALSE, space=0, freq=FALSE, main="Multiple histogram", xlab="x", ylab="Density")
	multhist(xnew, beside=TRUE, breaks=mh$breaks, freq=FALSE, main="Multiple histogram", xlab="x", ylab="Density")
	hist(xnew[[1]], breaks=mh$breaks, freq=FALSE, col=grey.colors(2)[1], main="Variable x, group 1", xlab="x1")
	hist(xnew[[2]], breaks=mh$breaks, freq=FALSE, col=grey.colors(2)[2], main="Variable x, group 2", xlab="x2")
par(mfrow=c(1,1))

mh$breaks

# MASS

dev.new()
ldahist(x,g,sep=TRUE, col=grey.colors(2))

dev.new()
ldahist(x,g,sep=FALSE,type="both")

ldahist(x,g,sep=FALSE,breaks=seq(-2.5, 2.5, 0.25))

ldahist(x,g,sep=FALSE)
