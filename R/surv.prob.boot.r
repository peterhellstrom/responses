# Bootstrap
q.fun <- function(x, na.rm=TRUE) quantile(x, c(0.025,0.5,0.975), na.rm=na.rm) # Function for extracting quantiles
randomSample <- function(df, n=nrow(df), replace=TRUE) df[sample(nrow(df), n, replace=replace),]

surv.prob.boot <- function(x, nboot=2000, na.rm=TRUE) {
	
	# Resample data, generate matrix with index variables (columns are replicates, rows are broods)
	inds <- sapply(1:nboot, function(i) sample(1:nrow(x), replace=TRUE, size=nrow(x)))
	# Point estimation for resampled data
	est <- lapply(1:nboot, function(i) surv.prob(x = x[inds[,i],]))
	# Extract data
	# Extract survival estimates (S.t) and survival function (S.f)
	S.t <- t(sapply(est, "[[", "S.t"))
	S.f <- t(sapply(est, "[[", "S.f"))
	
	# Calculate quantiles
	if (dim(x)[2] > 2) {
		out.boot <- list(
			S.t = t(apply(S.t, 2, q.fun, na.rm=na.rm)), S.f = t(apply(S.f, 2, q.fun, na.rm=na.rm)),
			SE.S.t = apply(S.t, 2, sd, na.rm=na.rm), SE.S.f = apply(S.f, 2, sd, na.rm=na.rm),
			nboot = nboot)
	} else {
		out.boot <- list(
			S.t = t(apply(S.t, 1, q.fun, na.rm=na.rm)), S.f = t(apply(S.f, 1, q.fun, na.rm=na.rm)),
			SE.S.t = apply(S.t, 1, sd, na.rm=na.rm), SE.S.f = apply(S.f, 1, sd, na.rm=na.rm),
			nboot = nboot)
	}
	out.boot
}

surv.prob.boot.n <- function(x, g, nboot=2000, na.rm=TRUE) {
	x.g <- split(x, g)
	lapply(x.g, surv.prob.boot, nboot=nboot, na.rm=na.rm)
}
