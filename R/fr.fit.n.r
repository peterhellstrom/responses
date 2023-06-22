# Fit many functional response datasets
# Wrapper function for fr.fit

# Must be updated to include start values in fr.fit!!!
# This function DOES NOT WORK AT ALL yet due to changes in fr.fit
#' @export
fr.fit.n <- function(data,method,eq) {

	nsim <- length(data)

	# Repeat the curve-fitting "nsim" number of times set:
	# Store the output in a list called res:
	# Output: x,y, AICc and Akaike weights, stored in res - a huge list.
	res <- lapply(1:nsim, function(i) try(fr.fit(data[[i]],method,eq),silent=TRUE))

	# Some generated data sets can not be fitted, so keep track of those:
	# rows is the number of succesful simulations
	# scc.res is a subset of res, containing only the succesful iterations:
	rows <- length(res[sapply(res, function(x) !inherits(x, "try-error"))])
	scc.res <- res[sapply(res, function(x) !inherits(x, "try-error"))]

	n.success <- which(sapply(res, function(x) !inherits(x, "try-error"))==TRUE)
	n.fail <- which(sapply(res, function(x) !inherits(x, "try-error"))==FALSE)

	# Extract Akaike weights
	Weights <- t(sapply(1:rows, function(i) scc.res[[i]]$AkaikeWeights))

	# FINAL output
	out <- list(
		out=scc.res,
		AICc.weights = Weights,
		iterations = nsim,
		succesful.iterations = rows,
		n.success=n.success,
		n.fail=n.fail
	)
	return(out)
}
