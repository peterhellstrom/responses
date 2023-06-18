# PDF for three-parameter t-distribution (NOT truncated)
dt3 <- function(x,mu,sigma,nu) {
	lambda <- 1/sigma^2 # squared inverse scale parameter (precision)
	(gamma((nu+1)/2) / gamma(nu/2)) * ((lambda / (pi*nu))^(1/2)) * ((1 + ((lambda*(x-mu)^2) / nu))^(-(nu+1)/2))
}


# Generate random variates from a truncated three-parameter t-distribution
rt3tr <- function(n, spec="t", a, b, mu, sigma, nu) {
	if (spec != "t") stop("You can only use the t-distribution")
	# Standardize truncation points
	a.st <- (a - mu) / sigma
	b.st <- (b - mu) / sigma
	# Generate random variates
	x <- (rtrunc(n=n, spec="t", df=nu, a=a.st, b=b.st) * sigma) + mu
	x
}
