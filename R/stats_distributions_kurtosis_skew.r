# This function computes the estimate of kurtosis of the distribution
# of x, following Box 6.2 in Sokal and Rohlf
kurtosis <- function(x) {
	m <- mean(x)
	s <- sd(x)
	y <- x - m
	sumy4 <- sum(y^4)
	n <- length(x)
	( (n+1)*n/((n-1)*(n-2)*(n-3)) ) * sumy4/s^4 - 3*(n-1)^2/((n-2)*(n-3))
}

# This function computes the estimate of skewness of the distribution
# of x, following Box 6.2 in Sokal and Rohlf
skew <- function(x) {
	m <- mean(x)
	s <- sd(x)
	y <- x - m
	sumy3 <- sum(y^3)
	n <- length(x)
	( n/((n-1)*(n-2)) ) * sumy3/s^3
}
