extr.try <- function(x, fn) { 

	sets <- length(x)
	# rows is the number of succesful simulations
	# scc.res is a subset of res, containing only the succesful iterations:
	inds <- sapply(x, function(x) !inherits(x, "try-error"))
	scc.x <- x[inds]
	rows <- length(scc.x)
	n.success <- which(inds == TRUE)
	n.fail <- which(inds == FALSE)

	cat(paste("Length of input list: ", sets, "\n", sep=""))
	cat(paste("Number of successes: ", rows, "\n", sep=""))
	cat(paste("Failure: ", n.fail, "\n", sep=""))

	out <- sapply(scc.x, fn)
	out

}

