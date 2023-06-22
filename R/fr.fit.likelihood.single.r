# Likelihood functions for functional responses
# All functions have normally distributed residuals
# Thus, the functions minimizes sums-of-squares and should give very similar results to nls()

# UPDATE 2013-11-14
# TODO: Check more about sigma...
# Previous versions did not include sigma as an estimated parameter
# I got this idea from the course at Grimsö with Tom Hobbs in 2008,
# but since then Tom has changed his course materials to include sigma as an estimated parameter.

# The original sigma statements are still in the code, but out-commented

# The new changes in these functions also have effects on fr.fit.r

# Null model, intercept only
#' @export
LL.null <- function(x, a, sigma) {
	y.hat <- rep(mean(a),length(x))
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Type I, Linear response without intercept
#' @export
LL.l <- function(y, x, a, sigma) {
	y.hat <- (a*x)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Michaelis-Menten, type II, IIIa & type IIIb
# Nested models, fix parameter theta=1 for Type II, theta=2 for Type IIIa
#' @export
LL.mm <- function(y, x, a, b, theta, sigma) {
	#y.hat <- TypeIIIb(x, a, b, theta)
	y.hat <- fr.mm(x, a, b, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Using weights
#' @export
LL.mm.w <- function(y, x, weights, a, b, theta, sigma) {
	#y.hat <- TypeIIIb(x, a, b, theta)
	y.hat <- fr.mm(x, a, b, theta)
	# sigma <- sqrt(sum(weights*((y.hat-y)^2)) / length(x))
	return(-sum(dnorm(y,y.hat,sigma/sqrt(weights),log=TRUE)))
}

# Supply vector with parameters, not sure if I have ever used this approach
# Delete?
#' @export
LL.mm.p <- function(p, x, y, sigma) {
	#if (is.null(names(p))) {
		#a <- p[1]
		#b <- p[2]
		#theta <- p[3]
	#} else {
		#if (any(c("a","b","theta") %in% names(p)==FALSE)) {
			#stop("Parameter names do not match")
		#} else {
			#a <- p[["a"]]
			#b <- p[["b"]]
			#theta <- p[["theta"]]
		#}
	#}
	a <- p[1]
	b <- p[2]
	theta <- p[3]

	y.hat <- fr.mm(x, a, b, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}


# Likelihood-function for two groups, type IIIb function
# Input variable: x,y,g [indexing variable for group membership]
#' @export
LL.mm.2gr <- function(a1, a2, b1, b2, theta1, theta2, sigma) {
	g <- as.numeric(g)
	a <- c(a1,a2)[g]
	b <- c(b1,b2)[g]
	theta <- c(theta1, theta2)[g]
	y.hat <- TypeIIIb(x=x, a=a, b=b, theta=theta)
	# ssq <- sum((y.hat-y)^2)
	# n <- length(x)
	# sigma <- sqrt(ssq/n)
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Michaelis-Menten, type IVa & b
# Nested models, fix theta=1 for type IVa
#' @export
LL.mm4 <- function(y, x, a, b, d, theta, sigma) {
	# y.hat <- TypeIVb(x, a, b, d, theta)
	y.hat <- fr.iv(x, a, b, d, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Holling's parameterization, type II, type IIIa & type IIIb
# Nested models, fix parameter theta to 1 for Type II, 2 for Type IIIa
#' @export
LL.h <- function(y, x, a, h, theta, sigma) {
	# y.hat <- TypeIIIb.h(x, a, h, theta)
	y.hat <- fr.h(x, a, h, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}
