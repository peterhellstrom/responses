library(mvtnorm)  
library(emdbook)

mu <- c(0,0) # vector with means
sigma <- matrix(c(1,0.5,0.5,1),nrow=2) # covariance matrix

# Plot the multivariate normal distribution with curve3d from package emdbook.
curve3d(dmvnorm(c(x,y),mu,sigma), from=c(-3,-3), to=c(3,3), sys3d="persp", theta=30, phi=30)

curve3d(dmvnorm(c(x,y),mu,sigma), from=c(-3,-3), to=c(3,3), sys3d="contour")
