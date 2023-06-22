#' @export
plot.acc.curve <- function(obj, plot.type="normal", xlim=NULL, addPoints=TRUE) {

	x <- obj$Samples
	y <- obj$Genotypes
	A <- obj$ParameterEstimates

	if (is.null(xlim)) {
		# Make a plot showing a) the accumulation curve and b) the three different fitted models
		# Extended plot type means
		if (plot.type=="extended") { xlims <- c(0,A[1,1]*A[1,2]-A[1,2]); ylims <- c(0, A[1,1]) }
		else if (plot.type=="normal") { xlims <- c(0,max(x)); ylims <- c(0,max(y)) }
		else stop ("The plot type you specified is not valid!")
		}

	if (!is.null(xlim)) {
		# Make a plot showing a) the accumulation curve and b) the three different fitted models
		# Extended plot type means
		if (plot.type=="extended") { xlims <- xlim; ylims <- c(0, A[1,1]) }
		else if (plot.type=="normal") { xlims <- xlim; ylims <- c(0,max(y)) }
		else stop ("The plot type you specified is not valid!")
		}

		par(mfrow=c(1,2))

		plot(c(0,x), c(0,y), pch=16, bty="l", xlab="Sample number",
			ylab="Unique genotypes", type="l", lwd=2, lty=1,
			ylim=ylims, xlim=c(0,max(x)),
			main="Accumulation curve", font.lab=2, las=1)
			points(c(0,x), c(0,y), col="red", pch=16, cex=0.7)

		plot(x=c(0,xlims), c(0,ylims), pch=16, bty="l", xlab="Sample number",
			ylab="Unique genotypes", type="n", lwd=2, lty=1,
			ylim=ylims, xlim=xlims,
			main="Model estimates", font.lab=2, las=1)

			if (addPoints == TRUE) points(c(0,x), c(0,y), col="red", pch=16, cex=0.7)
			curve(A[1,1]*x/(A[1,2]+x), from=0, to=xlims[2], add=T, type="l", lty=1, lwd=2, col="black")
			curve(A[2,1]*(1-exp(A[2,2]*x)), add=T, type="l", lty=2, lwd=2, col="green")
			curve(A[3,1]-(A[3,1]*(1-(1/A[3,1]))^x), add=T, type="l", lty=3, lwd=2, col="blue")
			abline(v=length(x), lty=2) # Draw vertical line at number of samples

			legend("topleft",c("Kohn","Eggert","Chessel"), col=c("black","green","blue"), lwd=c(2,2,2), lty=c(1,2,3), bty="n", cex=0.8)

			par(mfrow=c(1,1))
}
