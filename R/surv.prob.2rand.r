# Observed data - Calculate survival
# Input columns with number of individuals alive at each time step and grouping variable g
# g can only contain two levels!
#' @export
surv.prob.2rand <- function(x, g, R=5000, plot=TRUE) {
	if (length(levels(g)) > 2) stop("It's only possible to compare two survival functions")

	# Test statistic (Pollock et al 1989)
	surv.0 <- surv.prob.n(x=x, g=g)
	S.t.0 <- t(sapply(surv.0, "[[", "S.t"))
	SE.t.0 <- t(sapply(surv.0, "[[", "SE.t"))
	# Calculate test-statistic
	# z-score for each individual conditional survival time
	z0.t <- diff(apply(S.t.0, 2, rev)) / colSums(SE.t.0)
	# z.t <- diff(S.t.0) / colSums(SE.t.0)
	d2.obs <- sum(z0.t^2)

	# Randomization procedure
	d2.rand <- numeric(R) # Store output

	for (i in 1:R) {
		# Resample the grouping vector:
		g.rand <- sample(g, replace=FALSE, size=length(g))
		# Calculate survival for randomized sample
		surv.rand <- surv.prob.n(x=x, g=g.rand)
		S.t.rand <- t(sapply(surv.rand, "[[", "S.t"))
		SE.t.rand <- t(sapply(surv.rand, "[[", "SE.t"))
		# Calculate test-statistic
		# z-score for each individual conditional survival time
		z.t <- diff(apply(S.t.rand, 2, rev)) / colSums(SE.t.rand)
		# z.t <- diff(S.t.rand) / colSums(SE.t.rand)
		d2 <- sum(z.t^2)
		d2.rand[i] <- d2
	}

	# Get p-value
	p.rand <- length(which(d2.rand >= d2.obs)) / R

	if (plot) {
		# Histogram of test statistic z
		hist(d2.rand)
		abline(v=d2.obs, lty=2, col=2, lwd=2)
	}

	#out <- list(d2.obs=d2.obs, d2.rand=d2.rand, p=p.rand)
	out <- list(d2.obs=d2.obs, p=p.rand)
	out
}
