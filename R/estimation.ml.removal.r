# REMOVAL METHODS for closed populations
# MAXIMUM-LIKELIHOOD ESTIMATION
# Removal method
# From Borchers et al 2002
# Code adapted from WiSP function point.est.rm
# Note that interval estimation is not included in this function!!!
# Input: vector of raw captures

ml.rm <- function(captures) {

	rs <- c(0,captures)[-length(c(0,captures))]
	ns <- captures
	nall <- sum(captures)
	n.occ <- length(ns)
	Es <- 1 # mean group size, set to 1

	# Functions
	llk <- function(x) {
		N <- exp(x[1]) + nall
		p <- exp(x[2])/(1 + exp(x[2]))
		Ns <- N - cumsum(rs)
		temp1 <- lgamma(Ns + 1) - lgamma(Ns - ns + 1) - lgamma(ns + 1)
		temp2 <- ns * log(p)
		temp3 <- (Ns - ns) * log(1 - p)
		llk <- -sum(temp1 + temp2 + temp3)
		return(llk)
	}
            
	transform.xtoNtheta <- function(x) {
		N <- exp(x[1]) + nall
		theta <- exp(x[2])/(1 + exp(x[2]))
		return(c(N, theta))
	}
            
	transform.Nthetatox <- function(x) {
		x <- c(log(x[1] - nall), log(x[2]/(1 - x[2])))
		return(x)
	}

	# MLE estimation
	startNtheta <- c(((1 + 1/n.occ) * nall), max((1/(1 + 1/n.occ)), .Machine$double.xmin))
	startx <- transform.Nthetatox(startNtheta)
	res <- nlm(llk, startx)
	Nthetahat <- transform.xtoNtheta(res$estimate)
	Nhat.grp <- round(Nthetahat[1])
	Nhat.ind <- round(Nhat.grp * Es)
	phat <- Nthetahat[2]
	log.Likelihood <- -res$minimum
	AIC <- -2 * log.Likelihood + 2 * length(phat)

	# Output
	pointest <- list(sample = captures, Nhat.grp = Nhat.grp, Nhat.ind = Nhat.ind, phat = phat, Es = Es, log.Likelihood = log.Likelihood, AIC = AIC, created = date())
	class(pointest) <- "point.est.rm"
	return(pointest)
}
