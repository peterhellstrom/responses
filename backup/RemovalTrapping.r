# REMOVAL TRAPPING INDEX
# Index for relative abundance
# Various methods compiled from the literature

# Terminology (compared to the ML functions):
# a2 = captures
# T = effort (in ml.ce)

# Input values can be single values, or a matrix/list(?) containing data for many
# trapping sessions.

# Requires four inputs, nontargets (a1), targets (a2), traps (T) & duration (D)
# T = number of traps
# D = duration of study, number of days (total uncorrected trap effort = T * D)
# a0 = number of traps that were not sprung
# a1 = sum of the number of sprung, but empty, traps and traps that caught a non-target species
# a2 = number of traps that caught an individual of the target species

# Index from Fitzgerald, Efford & Karl
# New Zealand Journal of Zoology, 31:167-184 2004
# If a1 = 0, the point estimate of Fitzgerald's index equals Caughley's index (see below)
# Estimation of SE assumes that captures per trap follow a Poisson distribution
# 2) Compare Fitzgerald's estimate to Caughley's (from his classic book, 1977)
# Caughley did not account for competing species, only the target species
# If a1 in Fitzgerald's index > 0, then Caughley's estimator is biased low
# NOTE: Fitzgerald's index is really equation 4 in Linn & Downton (1975)
# 3) Nelson & Clark's (1973) index (Journal of Mammalogy, 1973)
# This index takes number of checks/total trap time into account
# The formula below is a modification by Peter Hellstr?m
# New, additional parameter:
# cpt = checks per trap (number of times each trap was checked)
# 4) Uncorrected density index
# Standard approach in many studies


DI <- function(targets, nontargets=0, traps, duration, cpt=duration, rowNames=NULL){
	
		a1 <- nontargets
		a2 <- targets
		effort <- traps*duration
		correction <- duration*(a1 + a2)/(2*cpt)
		corrected.effort <- effort - correction
		a0 <- effort - a1 - a2
	
		index <- data.frame(
					fitzgerald = -100*(a2/(a1+a2))*log(a0/effort),
					caughley = -100*log(1-(a2/effort)),
					nelson = 100*a2/corrected.effort,
					uncorrected = 100*a2/effort
				)
	
		# Fitzgerald
		se.index <- data.frame(
						fitzgerald = 100*sqrt(((a1*a2*(log(a0/effort)^2))/(effort-a0)^3)+(a2^2/(a0*effort*(effort-a0))))
					)
		
		input <- data.frame(
					targets = targets,
					nontargets = nontargets,
					traps = traps,
					duration = duration,
					cpt = cpt,
					effort = effort,
					correction = correction,
					corrected.effort = corrected.effort)

		if (is.character(rowNames)==TRUE | is.numeric(rowNames)) {
			if (length(rowNames) != nrow(index)) stop("Rownames does not match")
			rownames(input) <- rownames(index) <- rownames(se.index) <- rownames
		}
		
		out <- list(
			input = input,
			index = index,
			se.index = se.index
			)
		
		out
	}

# Sources:
# Otis DL, Burnham KP, White GC, Anderson DR (1978) Statistical inference from capture data on closed populations. Wildl Monogr 62:1-135 (p. 44-)
# White GC, Anderson DR, Burnham KP, Otis DL (1982) Capture-recapture and removal methods for sampling closed populations. Los Alamos National Laboratory, Los Alamos, New Mexico

captures <- function(N, sessions, c.prob) {

	E.n <- numeric(sessions)
	capt.n <- numeric(sessions)
	not.caught <- numeric(sessions)
	pop.start <- numeric(sessions)
	n.start <- N
	
	if (length(c.prob) == 1) c.prob <- rep(c.prob, sessions)
	if (length(c.prob) != sessions) stop("c.prob must be a single probability, or a vector of the same length as the number of sessions")

	for (j in 1:sessions) {
		E.n[j] <- N * ((1 - c.prob[j])^(j-1)) * c.prob[j] # This is probably incorrect if c.prob varies between sessions?
		capt.n[j] <- rbinom(1, n.start, c.prob[j])
		pop.start[j] <- n.start
		not.caught[j] <- n.start - capt.n[j]
		n.start <- not.caught[j]
	}

	mat <- matrix(ncol=5,nrow=sessions)
	mat[,1] <- 1:sessions
	mat[,2] <- pop.start
	mat[,3] <- E.n
	mat[,4] <- capt.n
	mat[,5] <- not.caught
	colnames(mat) <- c("Occasion", "Pop Size At Start", "Exp Captures", "Captures", "Not yet caught")
	mat
}

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


