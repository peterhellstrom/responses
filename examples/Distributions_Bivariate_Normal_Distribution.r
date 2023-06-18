library(ellipse)
library(emdbook)
library(mnormt)

# The density of a bivariate normal distribution is given by:
dbinorm <- function(x, mu, sigma) {
	a <- t(x-mu) %*% solve(sigma) %*% (x-mu)
	exp(-a/2)/(sqrt(det(sigma))*2*pi)
}

# A similar command in mnormt
# dmnorm(x, mean=mu, varcov=sigma)

# Determinant of sigma = sigma[1,1]*sigma[2,2]-sigma[1,2]*sigma[2,1]

# mu = vector with means, sigma = covariance matrix
mu <- c(0, 0)
sigma <- matrix(c(sqrt(2), 1, 1, sqrt(2)), ncol=2)
x <- matrix(runif(10,-3,3), ncol = 2)

# My function dbinorm can not take vectors, so we have to use apply:
sapply(1:nrow(x), function(i) dbinorm(x = x[i,], mu = mu, sigma = sigma))
# Identical with
dmnorm(x = x, mean = mu, varcov = sigma)

# Use Bolker's curve3d in emdbook to draw the bivariate distribution.
# Change formula a little bit... (allow x & y inputs)
dbinorm <- function(x, y, mu, sigma) {
	z <- c(x, y)
	a <- t(z - mu) %*% solve(sigma) %*% (z - mu)
	exp(-a / 2) / (sqrt(det(sigma)) * 2 * pi)
}

# Contour
curve3d(dbinorm(x, y, mu, sigma), from = c(-3, -3), to = c(3, 3), n = c(101, 101), sys3d = "contour")
# Perspective
curve3d(dbinorm(x, y, mu, sigma), from = c(-3, -3), to = c(3, 3), sys3d = "persp", theta = -30, phi = 40)

# Draw random variates from a multinormal distribution
n <- 200
d <- rmnorm(n, mu, sigma)
# Draw ellipses (from package ellipse)
dev.new()
plot(d[,1], d[,2])
lines(ellipse(sigma, centre = mu, level = 0.50), lty = 2)
lines(ellipse(sigma, centre = mu, level = 0.90), col = 2, lty = 2, lwd = 2)

