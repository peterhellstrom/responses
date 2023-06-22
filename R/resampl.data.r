# This function randomizes the input data and runs the analysis a specified number of times.
# It uses the acc.curve() function

#' @export
resampl.data <- function(samples, genotypes, resamplings, replace=FALSE, Log=FALSE) {

	nS <- length(samples)
	nUG <- length(unique(genotypes))
	nR <- resamplings
	x <- 1:nS
	# Create a matrix with datasets
	z <- sapply(1:resamplings, function(i) sample(genotypes, nS, replace=replace, prob=NULL) )
	# Calculate the cumulative number of unique genotypes in the sample z
	y <- sapply(1:resamplings, function(i) cumunique(z[,i]))
	if (Log == TRUE) y <- log(y)
	rownames(y) <- x
	colnames(y) <- paste("Smp", 1:resamplings, sep="")

	yLong <- data.frame(
				genotypes = c(y),
				samples = rep(x, times=resamplings),
				resample = as.factor(rep(colnames(y), each=nS)))
	yGr <- groupedData(genotypes ~ samples|resample, data=yLong)

	out <- list(
		nS = nS,
		nUG = nUG,
		resamplings = nR,
		Log = Log,
		x = x,
		z = z,
		y = y,
		yGr = yGr)
	out
}
