# NOTE! This function does not account for irregular census intervals
# E.g. time series must have only a single deltat / frequency!
#' @export
gr <- function(N, method="r", plot=TRUE) {

	if (class(N) != "ts") stop("Object N must be a time series object")

	na.inds <- which(is.na(N))
	if (length(na.inds) > 0) message("Object N contains missing/NA values")

	g <- switch(method,
		lambda = N[-1]/N[-length(N)],
		r = diff(log(N))
	)

	# Time component of the variable g is wrong (it is associated with t instead of t-1)
	# and must therefore be changed
	g <- ts(g, start=start(N)[1], deltat=deltat(N))

	title.text <- switch(method,
		lambda = expression(bold(paste("Growth rate =", lambda))),
		r = expression(bold(paste("Growth rate = log(", lambda, ")")))
	)

	if (plot) {

		# LOESS for N
		loess.N <- loess(N ~ time(N))
		# LOESS for g
		loess.g <- loess(g ~ time(g))
		# LOESS for dd in g
		dd.N <- N[-length(N)]
		loess.dd <- loess(g ~ dd.N)
		pred.dd.N <- seq(min(N,na.rm=T),max(N,na.rm=T),length=101)
		pred.loess.dd <- predict(loess.dd, newdata=pred.dd.N)

		# Density dependence in growth rates
		lm.dd <- lm(g ~ N[-length(N)])
		coefs <- summary(lm.dd)$coefficients[,1]
		p <- summary(lm.dd)$coefficients[2,4]
		lm.text <- paste(round(coefs[1],6), " + ", round(coefs[2],6), "x, p =", round(p,6), sep="")

		par(mfrow=c(2,2))
		plot(N, las=1, type="n", main="Population counts", xlab="Time", ylab="Counts")
			lines(N, lty=2)
			points(N, pch=16, cex=1)
			points(time(N), predict(loess.N), type="l", col=2, lwd=2)
		plot(g,las=1,type="n",main=title.text,xlab="Time",ylab="Growth rate")
			lines(g,lty=2)
			points(g,pch=16,cex=1)
			points(time(g), predict(loess.g), type="l", col=2, lwd=2)
		plot(N[-length(N)],las=1,g,pch=16,cex=1,xlab="Population density",ylab="Growth rate",main="Density-dependence in growth rates")
			points(pred.dd.N,pred.loess.dd,col=2,lwd=2,type="l")
		acf(g, main="Auto-correlation in growth rates", na.action=na.pass)
		if (length(na.inds) > 0) title(sub="N contains NA values")
		par(mfrow=c(1,1))
	}
	# Print output
	invisible(list(TimeSeries=N, GrowthRate=g))
}
