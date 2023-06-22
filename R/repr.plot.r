# Plot reproduction data (yearly averages), based on data summarized with the pop.samp function
#' @export
repr.plot <- function(x, y, g, error=NULL, n, xlims=range(x, na.rm=T), ylims=NULL,
	xlab="Rodent density index", ylab="Population average", main="",
	mar = c(5,5.5,1,2), cex=1.3, line=3.5,
	points=TRUE, pch = as.numeric(g), legend=TRUE, ...) {

	if (is.null(ylims)) {
		if (is.null(error)) {
			ylims <- c(0,ceiling(max(y, na.rm=T)))
		} else {
			ylims <- c(0,ceiling(max(y + error, na.rm=T)))
		}
	}

	op <- par(mar=mar)
	plot(x, y, type="n", xlab=xlab, ylab="", main=main, font.lab=2, las=1, cex.axis=cex, cex.lab=cex, ylim=ylims, xlim=xlims, bty="l")
	mtext(text=ylab, side=2, line=line, cex=cex, font=2)

	if (!is.null(error)) segments(x, y - error, x, y + error)
	if (points) points(x, y, pch=pch, cex=cex, ...)
	if (legend) repr.plot.legend()
	par(op)
}


repr.plot.points <- function(x, y, pch, error=NULL, cex=1.3, ...) {
	if (!is.null(error)) segments(x, y - error, x, y + error)
	points(x, y, pch=pch, cex=cex, ...)
}

repr.plot.legend <- function(pos="bottomright", cex=1.3) legend(pos, c("1970s","2000s"), pch=c(16, 21), pt.bg=c("black","white"), cex=cex, lty=c(1,2), bty="n")
