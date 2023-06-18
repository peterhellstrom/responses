
# Update code to only include the multi.disc with minimization of the Box-Draper criterion!
# Procedure similar to current LLmulti.disc2
# LLmulti.disc only estimates response to one prey species
# LLmulti.disc2 estimates parameters to both species and is the preferred approach.

LLmulti.disc <- function(y, x1, x2, a1, a2, h1, h2) {
inds <- which(y<0) # do not allow negative values
if (length(inds)>0) {
  y <- y[-inds]
  x1 <- x1[-inds]
  x2 <- x2[-inds]
  } 

y.hat <- multi.disc(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2)$f1
  ssq <- sum((y.hat-y)^2)
  n <- length(y)
  sigma <- sqrt(ssq/n)
  return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

#####################
 
LLmulti.disc2 <- function(y1, y2, x1, x2, a1, a2, h1, h2) {
# I can't recall where this code came from! It is likely correct, but... 
  # Observed response
  y <- matrix(cbind(y1,y2),ncol=2,byrow=F)
  # Predicted response
  y.hat <- multi.disc(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2)
  mu <- matrix(cbind(y.hat$f1,y.hat$f2),ncol=2,byrow=F)
  # Residual matrix
  z <- t(y - mu)
  n <- length(mu)
  # Variance-covariance matrix of residuals
  Sigma <- (z %*% t(z))/n
  eS <- eigen(Sigma, sym = TRUE, EISPACK = TRUE)
  ev <- eS$values
  logdetS <- try(determinant(Sigma, logarithm = TRUE)$modulus)
  iS <- try(solve(Sigma))
  ssq <- diag(t(z) %*% iS %*% z)
  loglik = -(n * (log(2 * pi)) + logdetS + ssq)/2
  sum(loglik)

return(-sum(dmvnorm(x=y, mu=mu, Sigma=Sigma, log=TRUE)))
}

LLmulti.disc2m <- function(y1, y2, x1, x2, a1, a2, h1, h2, m1, m2) {
# I can't recall where this code came from! It is likely correct, but... 
  # Observed response
  y <- matrix(cbind(y1,y2),ncol=2,byrow=F)
  # Predicted response
  y.hat <- multi.disc.m(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2, m1=m1, m2=m2)
  mu <- matrix(cbind(y.hat$f1,y.hat$f2),ncol=2,byrow=F)
  # Residual matrix
  z <- t(y - mu)
  n <- length(mu)
  # Variance-covariance matrix of residuals
  Sigma <- (z %*% t(z))/n
  eS <- eigen(Sigma, sym = TRUE, EISPACK = TRUE)
  ev <- eS$values
  logdetS <- try(determinant(Sigma, logarithm = TRUE)$modulus)
  iS <- try(solve(Sigma))
  ssq <- diag(t(z) %*% iS %*% z)
  loglik = -(n * (log(2 * pi)) + logdetS + ssq)/2
  sum(loglik)

return(-sum(dmvnorm(x=y, mu=mu, Sigma=Sigma, log=TRUE)))
}

LLmulti.disc2bd <- function(y1, y2, x1, x2, a1, a2, h1, h2) {
  # Observed responses
  y <- matrix(c(y1,y2),ncol=2)
  # Predicted responses
  y.hat <- multi.disc(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2)
  mu <- matrix(c(y.hat$f1,y.hat$f2),ncol=2)
  z <- y - mu
  # Box-Draper criterion (see chapter 4 in Bates & Watts 1988).
  # Minimizes the determinant of the crossproduct of the residual matrix.
  # Does not return a likelihood...but works still
  return(prod(svd(y - mu, nu=0, nv=0)$d)^2)
}

LLmulti.joly <- function(y, x1, x2, a.tot, h1, h2, alpha) {

y.hat <- multi.joly(x1=x1, x2=x2, a.tot=a.tot, h1=h1, h2=h2, alpha=alpha)
# Remove NA's
inds <- which(is.na(y.hat))

if (length(inds)>0) {
y <- y[-inds]
x1 <- x1[-inds]
x2 <- x2[-inds]
y.hat <- y.hat[-inds] }

ssq <- sum((y.hat-y)^2)
n <- length(y.hat)
sigma <- sqrt(ssq/n)
return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}
