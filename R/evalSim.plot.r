# Compare actual observed data with another (simulated) data set.
# Draw histogram of simulated data, and add min, max, density of observed data.
#' @export
evalSim.plot <- function(x, x.sim) {

	x.eval.max <- 100 * (table(x.sim >= max(x)) / length(x.sim)) # % of simulated values larger than max(observed)
	x.eval.min <- 100 * (table(x.sim <= min(x)) / length(x.sim)) # % of simulated values smaller than min(observed)

	hist(x.sim, breaks=30, col="steelblue", xlab="x", ylab="Density", main="Simulation vs. data", freq=FALSE)
		lines(density(x.sim), lwd=2)
		lines(density(x),col=3, lwd=2, lty=1)
		abline(v=min(x), col=4, lwd=2, lty=2)
		abline(v=max(x), col=2, lwd=2, lty=2)
		legend("topleft", c("Max","Min","Density"), title="Observations", col=c(2,4,3), lwd=2, lty=c(2,2,1), bty="n", cex=0.75)
	title(sub=paste(x.eval.min[2], "&", as.numeric(x.eval.max[2]),"% of simulated random variates are outside observed range"))
}
