# Convert to repeated measures input:
# Input is a data frame or matrix where rows are broods and time are columns
#' @export
rep.bn <- function(x) {
	if (length(colnames(x)) == 0) stop("Please supply column names for input matrix")
	brood <- factor(1:nrow(x))
	rownames(x) <- brood
	z <- ncol(x)
	temp <- lapply(2:z, function(i) data.frame(Alive=x[,i], Dead=x[,i-1]-x[,i], brood=brood, time=colnames(x)[i]))
	temp <- do.call("rbind", temp)
	out <- list(
		surv = cbind(Alive=temp$Alive, Dead=temp$Dead),
		brood = factor(temp$brood),
		timestep = factor(temp$time)
	)
	out
}

#inp <- rep.bn(x)
