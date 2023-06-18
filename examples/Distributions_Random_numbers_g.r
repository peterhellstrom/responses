##################################################
# Method for generating random variates for a sample, with a categorical grouping variable

xm <- c(1,10,20) # vector of means
xs <- c(2,2,2) # vectors of sd's

# Generate sites
catgs <- 1:3 # categories of grouping variable
n <- 100000 # number of random variates to sample

g <- sample(catgs,n,replace=TRUE) # create a grouping variable
table(g)/length(g) # check grouping variable

# Sample with unequal probabilities, adding arguments prob
g <- sample(catgs,n,replace=TRUE, prob=c(1/4, 1/4, 1/2)) # create a grouping variable
table(g)/length(g) # check grouping variable

# Generate random(normal) variates, by subsetting
.x <- rnorm(n, xm[g], xs[g])
length(.x)

# Check that the input parameters are retrieved
tapply(.x,g,mean)
tapply(.x,g,sd)

# Histogram of joint distribution
dev.new(width=7,height=5)
hist(.x, breaks=60, col="lightgrey", freq=FALSE)
lines(density(.x), lwd=2, col=2)

# Histogram 
library(lattice)
histogram(~.x|g, breaks=30) # Lattice histogram

dev.new(width=12,height=4)
par(mfrow=c(1,3))
for(i in catgs) {
	hist(.x[g==i], breaks=30, col="lightgrey", freq=FALSE, xlab=paste("x",i,sep=""), main=paste("Distribution of x",i,sep=""))
	curve(dnorm(x=x,xm[i],xs[i]),col=2,add=T,lwd=2)
}
par(mfrow=c(1,1))


##################################################
# This approach is implemented in the custom function genRandom:
source("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions/DistributionsFunctions.r")

# Test the genRandom-function:

# Normal
.x <- genRandom(distr="norm", para=list(n=10000, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, print=TRUE)
tapply(.x$x,.x$g,mean)
tapply(.x$x,.x$g,sd)
hist(.x$x, col="lightgrey", breaks=60)

# Normal, with unequal sampling probabilities
.x <- genRandom(distr="norm", para=list(n=10000, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, prob=c(1/4,1/4,1/2), print=TRUE)
tapply(.x$x,.x$g,mean)
tapply(.x$x,.x$g,sd)
hist(.x$x, col="lightgrey", breaks=60, freq=FALSE)

# Use with pre-defined grouping variable
n <- 10000
catgs <- 1:3 
#catgs <- factor(letters[1:3]) # Categories must be integer OR factor!!!
g <- sample(catgs, n, replace=TRUE)

.x1 <- genRandom(distr="norm", para=list(n=n, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, g=g, prob=c(1/4,1/4,1/2), print=TRUE)
.x2 <- genRandom(distr="norm", para=list(n=n, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, g=g, prob=c(1/4,1/4,1/2), print=TRUE)

hist(.x1$x, col="lightgrey", breaks=60, freq=FALSE)
table(.x1$g==.x2$g)


.x11 <- genRandom(distr="norm", para=list(n=n, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, prob=c(1/4,1/4,1/2), print=TRUE)
.x22 <- genRandom(distr="norm", para=list(n=n, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, prob=c(1/4,1/4,1/2), print=TRUE)

table(.x11$g==.x22$g)


# Uniform
.x <- genRandom(distr="unif", para=list(n=10000, min=c(0,10,50), max=c(1,25,75)), catgs=3, print=TRUE)
tapply(.x$x,.x$g,mean)
tapply(.x$x,.x$g,sd)
hist(.x$x, col="lightgrey", breaks=30, freq=FALSE)

# Lognormal
.x <- genRandom(distr="lnorm", para=list(n=10000, meanlog=c(1,10,20), sdlog=c(2,2,2)), catgs=3, print=TRUE)
hist(log(.x$x), col="lightgrey", breaks=60, freq=FALSE)
tapply(log(.x$x),.x$g,mean)
tapply(log(.x$x),.x$g,sd)

# Binomial distribution
.x <- genRandom(distr="binom", para=list(n=10000, size=c(1,1), prob=c(0.1,0.95)), catgs=2, print=TRUE)
table(.x$x,.x$g)
table(.x$g)/length(.x$g)

# Poisson distribution
.x <- genRandom(distr="pois", para=list(n=10000, lambda=c(5,7,20)), catgs=3, print=TRUE)
tapply(.x$x,.x$g,mean)
tapply(.x$x,.x$g,var)
tapply(.x$x,.x$g,sd)
hist(.x$x, col="lightgrey", breaks=30, freq=FALSE)

# Exponential distribution
.x <- genRandom(distr="exp", para=list(n=10000, rate=c(1/10, 1/5, 1/2)), catgs=3, print=TRUE)
tapply(.x$x,.x$g,mean)
tapply(.x$x,.x$g,sd)
hist(.x$x, col="lightgrey", breaks=60, freq=FALSE)

##################################################
# Truncated distributions
# Syntax is different, calls to function rtrunc
# Name of distribution is not specified in distr argument (has to be "trunc"), but rather as an item in the para list, argument spec.
# If not, the function rtrunc will not work.

n <- 100000
catgs <- 1:3
g <- sample(catgs, n, replace=TRUE)
table(g) / n

# prob <- c(1/4,1/4,1/2)
prob <- rep(1/length(catgs), length(catgs))

.x <- genRandom(distr="trunc", para=list(n=n, spec="norm", a=c(-3,-2,-1), b=c(3,2,1), mean=rep(0.5,3), sd=rep(1,3)), 
	catgs=3, prob=prob, print=TRUE)

table(.x$g) / n

tapply(.x$x,.x$g,mean)
tapply(.x$x,.x$g,sd)

dev.new(width=12, height=4)
par(mfrow=c(1,3))
	hist(.x$x[.x$g==1], breaks=30, col="steelblue", freq=FALSE, main="Group 1", xlab="x1")
	curve(dtrunc(x,spec="norm",a=-3, b=3, mean=0.5, sd=1),add=T)
	abline(v=mean(.x$x[.x$g==1]), lty=2, col=2)
	abline(v=0.5, lty=2, lwd=2)
	
	hist(.x$x[.x$g==2], breaks=30, col="steelblue", freq=FALSE, main="Group 2", xlab="x2")
	curve(dtrunc(x,spec="norm",a=-2, b=2, mean=0.5, sd=1),add=T)
	abline(v=mean(.x$x[.x$g==2]), lty=2, col=2)
	abline(v=0.5, lty=2, lwd=2)
	
	hist(.x$x[.x$g==3], breaks=30, col="steelblue", freq=FALSE, main="Group 3", xlab="x3")
	curve(dtrunc(x,spec="norm",a=-1, b=1, mean=0.5, sd=1),add=T)
	abline(v=mean(.x$x[.x$g==3]), lty=2, col=2)
	abline(v=0.5, lty=2, lwd=2)
par(mfrow=c(1,1))

##################################################
# Truncated 3-parameter t-distribution

hist(rt3tr(n=10000, a=-2, b=3, mu=0, sigma=1, nu=10), breaks=30, xlim=c(-3,4), col="steelblue", main="Histogram of x", freq=FALSE)
gen.trun(c(-2,3), family="TF", type="both")
curve(dTFtr(x, mu=0, sigma=1, nu=10), from=-2, to=3, n=1001 ,add=T, lwd=2)
curve(dt3(x, mu=0, sigma=1, nu=10), from=-3, to=4, n=1001, add=T, col=3, lwd=2)
legend("topright", c("Truncated","Untruncated"), lwd=c(2,2), lty=c(1,1), col=c(1,3), bty="n", title="3 parameter t-distribution")

n <- 10000
catgs <- 1:3
g <- sample(catgs, n, replace=TRUE)
prob <- rep(1/length(catgs), length(catgs))

para <- list(n=n, spec="t", a=c(0,17,28), b=c(20,25,32), mu=c(10,20,30), sigma=rep(1,3), nu=rep(10,3))

.x <- genRandom(distr="t3tr", para=para, catgs=3, prob=prob, g=g, print=TRUE)

dev.new(width=12, height=4)
par(mfrow=c(1,3))
	hist(.x$x[.x$g==1], breaks=30, col="steelblue", freq=FALSE, main="Group 1", xlab="x1")
	hist(.x$x[.x$g==2], breaks=30, col="steelblue", freq=FALSE, main="Group 2", xlab="x2")
	hist(.x$x[.x$g==3], breaks=30, col="steelblue", freq=FALSE, main="Group 3", xlab="x3")
par(mfrow=c(1,1))

##################################################
