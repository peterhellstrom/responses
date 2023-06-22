# Function for simulating data sets:
# Generate prey density (x), calculate deterministic y, and error normally distributed
# error to y.det. Parameterization is Michaelis-Menten.

# Inputs:
# parms = functional response parameters a,b,c (must be a list)
# Functional response parameters:
# a = asymptote (miximum kill rate)
# b = half-saturation constant
# c = exponent controlling shape of curve
# lims = min and max x-values (must be called min.x and max.x, list)
# sigma = error (rse)
# n = number of data points in each set
# sets = generate "sets" number of data sets
# xmethod = how to generate x-values:
# full.x.range = draw n uniform values from min.x to max.x
# even.x.range = split x-range into 4 bins (determined by the quantiles of min.x and max.x)# , assign equal (or nearly equal if n isn't a multiple of 4) number of values in each bin.

#' @export
fr.sim.data <- function(parms,lims,sigma,n,sets,xmethod) {

out <- vector("list",sets)

for (i in 1:sets) {

if (xmethod=="full.x.range") {
  x <- runif(n=n, min=lims$min.x, max=lims$max.x)
  }

if (xmethod=="even.x.range") {

  groups <- quantile(c(lims$min.x, lims$max.x))
  x.vec <- rep(1:4,times=ceiling(n/4))
  n.samp <- table(sample(x.vec,n,replace=FALSE))

  x <- c(
        runif(n.samp[1],min=groups[1],max=groups[2]),
        runif(n.samp[2],min=groups[2],max=groups[3]),
        runif(n.samp[3],min=groups[3],max=groups[4]),
        runif(n.samp[4],min=groups[4],max=groups[5])
        )
  }

  y.hat <- parms$a*x^parms$c / (parms$b^parms$c + x^parms$c)
  # Draw random normal variates and take absolute value
  # Not the best solution, should find a better solution
  # based on Poission or negbin
  y <- abs(rnorm(x,mean=y.hat,sd=sigma))

  out[[i]] <- list(
  input=data.frame(a=parms$a, b=parms$b, c=parms$c, min.x=lims$min.x, max.x=lims$max.x,
        sigma=sigma, n=n, sets=sets),
  x=x,
  y.hat=y.hat,
  y=y)

  }
out
}
