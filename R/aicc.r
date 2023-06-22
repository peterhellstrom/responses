# Stand alone version for calculating AIC and AICc values
# Use this function to compare AICc values
# AICc-functions are also available in the packages bbmle, MuMIn, AICcmodavg

# NOTE: beware that the attribute nobs differs between packages and functions
# nls & gnls return the actual number of observations
# gls return the number of observations - the number of model paramters (not including variance(s))!

#' @export
aicc <- function(object, nobs=NULL) {
	# Get the number of datapoints and estimated parameters for different object classes
	# object.class <- class(object)[1]
	# cat(paste("class attribute of object:",object.class), "\n")
	if (!is.null(nobs)) {
		if(length(nobs) > 1) stop("nobs must be a scalar!")
		n <- nobs
	} else {
		n <- attr(logLik(object), "nobs")
	}
	k <- attr(logLik(object), "df")
	aic <- -2*logLik(object) + 2*k # AIC
	penalty <- (2*k*(k+1)) / (n-k-1) # penalty term
	if (n-k-1 == 0) cat("Warning: denominator of penalty term equals 0.\nAICc is infinite.\n")
	AICc <- aic + penalty
	# Create output vector
	out <- c('logLik'=logLik(object), 'k'=k, 'n'=n, 'AIC'=aic, 'AICc'=AICc, 'penalty'=penalty)
	out
}

# Extract AICc info and calculate statistics for many data sets
#' @export
aicc.n <- function(object, sort=TRUE, round=FALSE, nobs=NULL, digits=5) {
	# Input is a list with objects (i.e. fitted models)
	n <- length(object)

	if (!is.null(nobs)) {
		tab <- data.frame(t(sapply(object,aicc,nobs)))
		tab$n <- nobs
	} else {
		tab <- data.frame(t(sapply(object,aicc)))
	}
	rownames(tab) <- names(object)
	index <- which(tab[,"AICc"] == min(tab[,"AICc"])) # Find model with smallest AICc

	deltai <- tab[,"AICc"] - min(tab[,"AICc"]) # AICc-delta values
	rel.like <- exp(-deltai / 2) # Relative likelihood
	wi <- rel.like / sum(rel.like) # Akaike weights
	ER <- wi[index] / wi # Information ratio
	ranking <- rank(deltai)

	# Create output (data frame)
	out <- data.frame(
	  'logLik' = tab$logLik,
	  'K' = tab$k,
	  'n'=tab$n,
	  'AIC' = tab$AIC,
	  'AICc' = tab$AICc,
		'delta' = deltai,
		'weigths' = wi,
		'ER' = ER,
		'rank' = ranking)

	rownames(out) <- names(object)

	if (sort) out <- out[order(out[,"rank"]),]
	if (round) {
		round(out,digits)
	} else
	out
}
