lnorm2norm <- function(mu,sigma) {
	names(mu) <- names(sigma) <- NULL
	 m <- exp(mu + (sigma^2)/2)
	 v <- exp(2*mu + sigma^2)*(exp(sigma^2)-1)
	 sd <- sqrt(v)
	 c(m=m, sd=sd, v=v)
}

norm2lnorm <- function(m,sd) {
	names(m) <- names(sd) <- NULL
	v <- sd^2
	mu <- log((m^2) / sqrt(v + m^2))
	sigma <- sqrt(log(v/(m^2) + 1))
	c(mu=mu,sigma=sigma,v=v)
}
