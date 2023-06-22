#' @export
fr.sim.data.multi <- function(parms,lims,sigma,n,sets,cor.method,...) {

# Simulate two-prey systems, mainly for testing estimation functions
# Works like:
# 1) create prey densities, draw from uniform distribution, or from a multivariate normal
# (option cor.method). Drawing from a multivariate normal distribution allows for correlation
# between prey numbers, but a variance-covariance matrix has to be defined for use.
# 2) Calculate the deterministic (y.hat) value, based on the values in the parms argument
# 3) Calculate the stochastic value (observation error) based on y.hat ("mean"), by drawing from
# a normal distribution with sd=sigma argument.
# cova = covariance matrix???

out <- vector("list",sets)

for (i in 1:sets) {

  if (cor.method=="null") {

  x1 <- runif(n=n, min=lims$min.x1, max=lims$max.x1)
  x2 <- runif(n=n, min=lims$min.x2, max=lims$max.x2)

  }

  if (cor.method=="correlated") {

  cova <- cova

  x <- mvrnorm(n, mu = c(cor.parms$x1,cor.parms$x2), Sigma = cova)

  inds <- which(x[1:n,]<0)
  x[inds] <- NA

  x1 <- x[,1]
  x2 <- x[,2]
  }

  y.hats <- multi.disc(x1=x1,x2=x2,a1=parms$a1,a2=parms$a2,h1=parms$h1,h2=parms$h2)

  y1.hat <- y.hats$f1
  y1 <- abs(rnorm(x1,mean=y1.hat,sd=sigma))
  # inds1 <- which(y1<0)
  # y1[inds1] <- NA

  y2.hat <- y.hats$f2
  y2 <- abs(rnorm(n,mean=y2.hat,sd=sigma))
  # inds2 <- which(y2<0)
  # y2[inds2] <- NA

  out[[i]] <- list(
  input = data.frame(a1=parms$a1, a2=parms$a2, h1=parms$h1, h2=parms$h2, min.x1=lims$min.x1, max.x1=lims$max.x1,
          min.x2=lims$min.x2, max.x2=lims$max.x2, sigma=sigma, n=n, sets=sets),
  x1=x1,
  y1.hat=y1.hat,
  y1=y1,
  x2=x2,
  y2.hat=y2.hat,
  y2=y2)

  }
out
}
