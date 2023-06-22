# Plot a gamlss.tr object and the input data
# object = name of gamlss-object.
# a & b = lower resp. upper truncation point
# breaks = number of bins for histogram.
#' @export
plot.Gamlss <- function(object, a=NULL, b=NULL, breaks=50) {

	x <- as.numeric(object$y)

	if (is.null(a)) a <- min(x)
	if (is.null(b)) b <- max(x)

	n <- length(x)
	distr <- object$family[1]
	pars.names <- rownames(coef.Gamlss(object))
	pars <- as.numeric(coef.Gamlss(object)$transformed)
	names(pars) <- pars.names

	pars.str <- paste(pars.names,"=",pars)
	pars.str <- paste(rep(pars.str),collapse=", ")

	qstr <- paste("q",distr,"(p=ppoints(n),", pars.str, ")", sep="")
	qx <- eval(parse(text=qstr))
	qx.resid <- residuals(object) # ?residuals.gamlss, GAMLSS returns quantile residuals
	qx.resid <- qx.resid[!is.infinite(qx.resid)] # Remove infinite values

	dstr <- paste("d",distr,"(x,", pars.str, ")", sep="")
	dstr <- paste("curve(expr=", dstr, ", from=a, to=b, n=1001, add=T, col=2, lwd=2)" ,sep="")

	par(mfrow=c(2,2))
		hist(x, col="steelblue", breaks=breaks, freq=FALSE, main=paste("f:",distr))
		lines(density(x))
		eval(parse(text=dstr))

		plot(qx, sort(x), xlab="Theoretical quantiles", ylab="Sample quantiles", main=paste("qq:",distr), cex=0.75)
		abline(0,1,lty=2)

		hist(qx.resid, col="steelblue", breaks=breaks, freq=FALSE, xlab="Quantile residuals", main=paste("resid(",distr, ")", sep=""))
		curve(dnorm(x), n=1001, col=2, lwd=2, add=T)

		qqnorm(qx.resid, xlab="Theoretical quantiles", ylab="Sample quantiles", main=paste("qq: resid(",distr, ")", sep=""), cex=0.75)
		qqline(qx.resid, lty=2)
	par(mfrow=c(1,1))

}
