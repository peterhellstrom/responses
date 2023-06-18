################################################################################
# Inverse transform sampling, Inversion method
# see http://en.wikipedia.org/wiki/Inverse_transform_sampling
# This example from http://stats.stackexchange.com/questions/12843/generating-random-samples-from-a-custom-distribution
################################################################################

# Original function, pdf & cdf

f.fun <- function(x) 1-x^2
curve(f.fun(x),from=0,to=1)
integrate(f.fun,lower=0,upper=1) # scaling factor for pdf = 1/(2/3) = 3/2

# Convert to pdf
f.pdf <- function(x) (3/2)*(1-x^2)
integrate(f.pdf,lower=0,upper=1) # Must be 1!

# Integrate the pdf to get the cdf
f.cdf <- function(x) (3/2)*(x-(x^3)/3)

# Plot pdf & cdf
dev.new(width=12,height=6)
par(mfrow=c(1,2))
curve(f.pdf(x),main="pdf",xlim=c(0,1))
curve(f.cdf(x), main="cdf",xlim=c(0,1))
par(mfrow=c(1,1))

# Generate random numbers
nsamples <- 100000
x <- runif(nsamples) 
f <- function(x, u) (3/2)*(x-(x^3)/3) - u 

# Illustrate the function f (check various values of u in the range 0 < u < 1.
# It must be possible to find a root in the x-range 0 to 1 (in this case...).
curve(f(x=x,u=0.5),from=0,to=1)
abline(h=0, lty=2, col=2)

# Use uniroot to find the root values (the random values)
# Define the function, note the extra argument 'u = x' in the call to uniroot.
my.uniroot <- function(x) uniroot(f, c(0, 1), tol = 0.0001, u = x)$root
# Use vapply to "loop" over all uniform values
system.time(r <- vapply(x, my.uniroot, numeric(1)))

hist(r,col="steelblue",breaks=30,freq=FALSE) # plot simulated random numbers
curve(f.pdf(x),add=TRUE,col=2,lwd=2) # add pdf to check that simulation was successful

################################################################################
# Some other examples (AR(2)-space)

f.fun <- function(x) 1-(x^2)/4
curve(f.fun(x),from=0,to=2)
integrate(f.fun,lower=0,upper=2) # scaling factor = 1/(4/3) = 3/4

f.pdf <- function(x) (3/4)*(1-x^2/4)
integrate(f.pdf,lower=0,upper=2) # must be = 1

f.cdf <- function(x) -(x*(-12 + x^2))/16

dev.new(width=12,height=6)
par(mfrow=c(1,2))
curve(f.pdf(x),main="pdf",xlim=c(0,2))
curve(f.cdf(x), main="cdf",xlim=c(0,2))
par(mfrow=c(1,1))

# Generate random numbers
nsamples <- 100000
x <- runif(nsamples) 
f <- function(x, u) (-(x*(-12 + x^2))/16) - u 

my.uniroot <- function(x) uniroot(f, c(0, 2), tol = 0.0001, u = x)$root 
system.time(r1 <- vapply(x, my.uniroot, numeric(1)))

hist(r1,col="steelblue",breaks=30,freq=FALSE)
curve(f.pdf(x),from=0,to=2,add=TRUE,col=2)
################################################################################

f.fun <- function(x) 2*sqrt(1-x)
curve(f.fun(x),from=0,to=1)
integrate(f.fun,lower=0,upper=1) # scaling factor = 1/(4/3) = 3/4

f.pdf <- function(x) (3/4)*(2*sqrt(1-x))
integrate(f.pdf,lower=0,upper=1) # must be = 1

f.cdf <- function(x) sqrt(1-x)*(x-1)+1

dev.new(width=12,height=6)
par(mfrow=c(1,2))
curve(f.pdf(x),main="pdf",xlim=c(0,1))
curve(f.cdf(x), main="cdf",xlim=c(0,1))
par(mfrow=c(1,1))


nsamples <- 100000
x <- runif(nsamples) 
f <- function(x, u) (sqrt(1-x)*(x-1)+1) - u 

my.uniroot <- function(x) uniroot(f, c(0, 1), tol = 0.0001, u = x)$root 
system.time(r2 <- vapply(x, my.uniroot, numeric(1)))

hist(r2,col="steelblue",breaks=30,freq=FALSE)
curve(f.pdf(x),from=0,to=1,add=TRUE,col=2)

ar2.plot(k=NULL,v=NULL)
points(r1,r2-1,pch=".",col=2)

################################################################################
