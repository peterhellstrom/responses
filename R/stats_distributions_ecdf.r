# Plot cumulative density function
# Possible to plot several distributions in the same plot, by using a list:

ecdfn <- function(x) {

if (class(x) == "matrix") {
	sortx <- sapply(1:ncol(x), function(i) sort(x[,i]))
	ecdfx <- sapply(1:ncol(x), function(i) (1:length(x[,i])) / length(x[,i]))
	
	
	plot(sortx, ecdfx, type="n", xlab="x", ylab="Fn(x)")
	abline(h=c(0,1), lty=2, col="grey")
	if (!is.null(colnames(x))) legend("bottomright", colnames(x), pch=16, cex=1, col=1:length(x), bty="n")
	points(sortx, ecdfx, pch=16, cex=0.7, col=rep(1:ncol(x),each=nrow(x)))
	for (i in 1:ncol(x)) lines(sortx[,i], ecdfx[,i], col=i)
}

if (class(x) == "list") {
	sortx <- lapply(1:length(x), function(i) sort(x[[i]]))
	ecdfx <- lapply(1:length(x), function(i) (1:length(x[[i]])) / length(x[[i]]))
	
	
	plot(unlist(sortx), unlist(ecdfx), type="n", xlab="x", ylab="Fn(x)")
	abline(h=c(0,1), lty=2, col="grey")
	if (!is.null(names(x))) legend("bottomright", names(x), pch=16, cex=1, col=1:length(x), bty="n")
	for (i in 1:length(x)) points(sortx[[i]], ecdfx[[i]], pch=16, cex=0.7, col=i)
	for (i in 1:length(x)) lines(sortx[[i]], ecdfx[[i]], col=i)
	}

	out <- list(x=x, sortx=sortx, ecdfx=ecdfx)
	out
}

# Another function for plotting an ecdf
ecdf2 <- function(x, plot=TRUE) {
	sortx <- sort(x)
	ecdfx <- (1:length(x)) / length(x)
	if (plot == TRUE) {
		
			plot(sortx, ecdfx, type="n", xlab="x", ylab="Fn(x)")
			abline(h=c(0,1), lty=2, col="grey")
			lines(sortx, ecdfx)
			points(sortx, ecdfx, pch=16, cex=.7)
		}
	out <- list(x=x,sortx=sortx,ecdfx=ecdfx)
	out
}
