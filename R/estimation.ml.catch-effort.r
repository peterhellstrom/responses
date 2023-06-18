# REMOVAL METHODS for closed populations
# MAXIMUM-LIKELIHOOD ESTIMATION
# Catch-effort method

# NOT WORKING GOOD FOR LARGE NUMBERS!!!!
# Solution: rescale effort
# From Borchers et al 2002
# Code adapted from WiSP function point.est.ce
# Note that interval estimation is not included in this function!!!
# Input: vectors containing raw captures and effort per trapping occasion

ml.ce <- function(captures, effort, iterlim=1000) {

	rs <- c(0,captures)[-length(c(0,captures))]
	ns <- captures
	nall <- sum(captures)
	n.occ <- length(ns)
	Es <- 1 #mean group size, set to 1
	effort <- effort

	# Functions
	llk <- function(x) {
		N <- exp(x[1]) + nall
		theta <- exp(x[2])
		Ns <- N - cumsum(rs)
		ps <- 1 - exp(-theta * effort)
		temp1 <- lgamma(Ns + 1) - lgamma(Ns - ns + 1) - lgamma(ns + 1)
		temp2 <- ns * log(ps)
		temp3 <- (Ns - ns) * log(1 - ps)
		llk <- -sum(temp1 + temp2 + temp3)
		return(llk)
	}

	transform.xtoNtheta <- function(x) {
		N <- exp(x[1]) + nall
		theta <- exp(x[2])
		return(c(N, theta))
	}

	transform.Nthetatox <- function(x) {
		x <- c(log(x[1] - nall), log(x[2]))
		return(x)
	}

	# MLE estimation
	startNtheta <- c(((1 + 1/n.occ) * nall), (ns[1]/nall))
	startx <- transform.Nthetatox(startNtheta)
	res <- nlm(llk, startx, iterlim=iterlim, gradtol = 1e-15)
	Nthetahat <- transform.xtoNtheta(res$estimate)
	Nhat.grp <- Nthetahat[1]
	theta <- Nthetahat[2]
	log.Likelihood <- -res$minimum
	AIC <- -2 * log.Likelihood + 2 * length(theta)
	pshat <- 1 - exp(-theta * effort)
	Nhat.grp <- round(Nhat.grp)
	Nhat.ind <- round(Nhat.grp * Es)

	# Output
	pointest <- list(sample = captures, effort = effort, Nhat.grp = Nhat.grp, Nhat.ind = (Nhat.grp * Es), theta = theta, phat = pshat, Es = Es, log.Likelihood = log.Likelihood, 
	AIC = AIC, created = date())
	class(pointest) <- "point.est.ce"
	return(pointest)
}
