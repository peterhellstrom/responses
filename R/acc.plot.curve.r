# Plot the population size estimates (min,max & mean) from acc.curve.resamp(),
# and also the distribution of the population size estimate 'a'.
# Requires a stored object fitted with acc.curve.resamp().
#' @export
plot.curve <- function(obj, model, ylim=NULL) {

	ests <- obj$Coefficients[,1:5]

	ests.mean <- apply(ests, 2, mean)
	ests.median <- apply(ests, 2, median)

	# ests.min <- apply(ests, 2, min)
	# ests.max <- apply(ests, 2, max)
	ests.min <- sapply(1:ncol(ests), function(i) quantile(ests[,i],c(0.025)))
	ests.max <- sapply(1:ncol(ests), function(i) quantile(ests[,i],c(0.975)))

	if (model == "Kohn") inds <- c(1,4)
	if (model == "Eggert") inds <- c(2,5)
	if (model == "Chessel") inds <- 3

	xlims <- c(0, obj$SampleSize)
	if (is.null(ylim)) ylims <- c(0, ests.max[inds][1])
	if(!is.null(ylim)) ylims <- ylim

	plot(x=xlims, y=ylims, type="n", xlab="Sample number", ylab="Discovered genotypes",
		font.lab=2, las=1, bty="l", main=paste("Population size based on ", model, " equation", sep=""))

	abline(h=ests.min[inds][1], lty=2, col=2)
	abline(h=ests.max[inds][1], lty=2, col=2)
	abline(h=ests.median[inds][1], lty=1, col=2)

	if (model != "Chessel") {
		str.min <- parse(file="", text=paste("curve(",model, "(x=x, a=ests.min[inds][1], b=ests.min[inds][2]), n=101, lty=2, add=T)", sep=""))
		str.max <- parse(file="", text=paste("curve(",model, "(x=x, a=ests.max[inds][1], b=ests.max[inds][2]), n=101, lty=2, add=T)", sep=""))
		str.median <- parse(file="", text=paste("curve(",model, "(x=x, a=ests.median[inds][1], b=ests.median[inds][2]), n=101, lty=1, add=T)", sep=""))

		eval(str.min)
		eval(str.max)
		eval(str.median)
		}

	if (model == "Chessel") {
		curve(Chessel(x, a=ests.min[inds]), from=0, to=obj$SampleSize, n=101, lty=2, add=T)
		curve(Chessel(x, a=ests.max[inds]), from=0, to=obj$SampleSize, n=101, lty=2, add=T)
		curve(Chessel(x, a=ests.median[inds]), from=0, to=obj$SampleSize, n=101, lty=1, add=T)
	}

}
