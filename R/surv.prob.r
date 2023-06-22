# surv.prob
# KAPLAN-MEIER product-type estimator
# See also flint.km.r for revised versions of some of the functions
# After Flint et al. 1995
# TO DO:
# Add more error handling

# Kaplan-Meier type: Calculations for standard errors: Should censored broods be included or not???

# Known errors in Kaplan-Meier type function:
# 1) Code breaks if all broods in a bootstrap sample are censored (should be uncommon, but happens for small samples)
# Some errors were caused by me trying to name columns in bootstrap output.
# Now bootstrap output are given without column names corresponding to time.
# When bootstrapping with few observations per time unit a common warning is: In qt(p, df, lower.tail, log.p) : NaNs produced

# GENERAL POINTS OF IMPORTANCE:
# - Non-independence among siblings
# - Hatching order (possibly) influences survival probability
# - Non-constant probability of survival over time, i.e. if one chick dies, the probability of survival for the remaining chicks may increase or decrease.
# - Success depends on number of "trials"?
# - Consider a death as a "removal" - will survival probabilities follow a hypergeometric or multinomial distribution?
# - beta-binomial distribution for survival probability?, see Bolker's book!

#' @export
surv.prob <- function(x) {

	# x is a matrix with broods as rows and time as columns
	# z = number of time steps, S = survival (proportion), w = weights, numbers at risk
	# censored, non.censored
	# S.t = weighted mean of the estimates for individual broods, S.f = Kaplan-Meier survival function estimates of the proportion surviving from day 0,
	# SE.it, SE.t,
	# n.bar = Mean brood size at time t, M = number of broods (or marked females with broods)

	z <- ncol(x) # Number of time steps = z-1
	counts <- nrow(x) # Number of rows/broods

	names.list <- list('Brood' = rownames(x), 'Time' = paste(colnames(x)[1:z-1], colnames(x)[2:z], sep="-"))

	# Calculate survival
	S <- sapply(1:(z-1), function(i) x[,i+1] / x[,i]); dimnames(S) <- names.list
	# Numbers at risk, used as weights.
	w <- cbind(x[,1:(z-1)]); dimnames(w) <- names.list

	# Logical, look for undefined survival estimates in each brood
	for (i in 2:z) undefined <- length(which(is.nan(S)))

		if (undefined > 0) {
			# Remove undefined values, using the any(), is.nan(), and drop functions
			x <- as.matrix(x[apply(S*w, 1, function(x)!any(is.nan(x))),, drop=F])
			w <- as.matrix(w[apply(S*w, 1, function(x)!any(is.nan(x))),, drop=F])
			S <- as.matrix(S[apply(S, 1, function(x)!any(is.nan(x))),, drop=F])
		}

		# CHANGE HERE!
		# If a brood only has a single observation (NA), survival can't be calculated and that brood should be removed:
		# nas <- apply(S, 1, function(x) length(which(is.na(x) & !is.nan(x))))
		nas <- apply(S, 1, function(x) length(which(is.na(x)))) # Count number of NA's on each row
		nas.inds <- which(nas == ncol(S))

		if (length(nas.inds > 0)) {
			x <- as.matrix(x[-nas.inds,])
			S <- as.matrix(S[-nas.inds,])
			w <- as.matrix(w[-nas.inds,])
		}

		valid.counts <- dim(S)[1] - sapply(1:ncol(S), function(i) length(which(is.na(S)[,i])))
		names(valid.counts) <- names.list$Time
		censored <- which(is.na(S*w))
		c.count <- length(censored)

		if (c.count == 0) {
			S.t <- if (z == 2) (sum(S*w, na.rm=T) / sum(w)) else (colSums(S*w, na.rm=T) / colSums(w, na.rm=T))
			S.f <- cumprod(S.t)
			SE.it <- w^2 * t((t(S)-S.t)^2)
			n.bar <- colMeans(as.matrix(x[,-c(z)]), na.rm=T)
			M <- valid.counts
			SE.denom <- M * n.bar[-z]^2 * (M-1)
			SE.t <- sqrt(colSums(SE.it, na.rm=T) / SE.denom[1:(z-1)])
			non.censored <- valid.counts
		} else if (c.count > 0) {
			if (z == 2) {
				S.t <- sum(S*w, na.rm=T) / sum(w[-censored])
				S.f <- cumprod(S.t)
				SE.it <- w^2 * t((t(S)-S.t)^2)
				# CODE STOPS HERE if bootstrap sample ends up with all samples being censored.
				n.bar <- colMeans(x[-censored,-z], na.rm=T)
				M <- valid.counts - length(censored)
				SE.denom <- M * n.bar[-z]^2 * (M-1)
				SE.t <- sqrt(colSums(SE.it, na.rm=T) / SE.denom[1:(z-1)])
				non.censored <- valid.counts - length(censored)
			} else if (z > 2) {
				w <- (S*w)/S
				S.t <- colSums(S*w, na.rm=T) / colSums(w, na.rm=T)
				S.f <- cumprod(S.t)
				SE.it <- w^2 * t((t(S)-S.t)^2)
				xS <- ifelse(!is.na(S),1,NA) # A matrix with zero's and ones indicating valid, non-censored observations (==1)
				M <- apply(xS,2,sum,na.rm=T)
				n.bar <- colMeans(xS * x[,-z], na.rm=T)
				SE.denom <- M * n.bar^2 * (M-1)
				SE.t <- sqrt(colSums(SE.it, na.rm=T) / SE.denom[1:(z-1)])
				non.censored <- length(which(!is.na(apply(xS,1,sum))))
			}
		}

		# Confidence limits
		cil.S.t <- S.t - (qt(0.975, df=M) * SE.t)
		ciu.S.t <- S.t + (qt(0.975, df=M) * SE.t)

		# output list
		out <- list(
			x = x, S=S, w=w,
			S.t = S.t, S.f = S.f, SE.it = SE.it, SE.t = SE.t, CI = rbind(CIL = cil.S.t, CIU = ciu.S.t),
			n.bar = n.bar, broods = counts, valid.broods = valid.counts,
			M = M, non.censored = non.censored, undefined = undefined)

		class(out) <- "flint"
		out
}

# x is a data frame with number of individuals in each brood
# g is a grouping variable (factor)
#' @export
surv.prob.n <- function(x, g) {
	x.split <- split(x, g)
	lapply(x.split, surv.prob)
}

# Calculate test statistic D^2 for difference between two samples [survival functions]
#' @export
d2.obs <- function(x, g) {
	if (length(levels(g)) > 2) stop("Only two groups can be compared!")
	surv.0 <- surv.prob.n(x=x, g=g)
	S.t.0 <- t(sapply(surv.0, "[[", "S.t"))
	SE.t.0 <- t(sapply(surv.0, "[[", "SE.t"))
	# Calculate test-statistic
	# z-score for each individual conditional survival time
	z0.t <- diff(apply(S.t.0, 2, rev)) / colSums(SE.t.0)
	# z.t <- diff(S.t.0) / colSums(SE.t.0)
	d2.obs <- sum(z0.t^2)
	list('S.t'=S.t.0, 'SE.t'=SE.t.0, 'z'=z0.t, 'd2'=d2.obs)
}

