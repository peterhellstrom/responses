# Make a plot that compares means and q-Intervals for the three different models
# Requires a stored object fitted with acc.curve.resamp().

plot.est.comp <- function(obj) {

	A <- obj$ParameterSummary

	dev.new(width=8,height=6)
	# Set plotting region and plot mean values
	plot(A["Mean",1:3], xlim=c(0.75, 3.25), ylim=c(min(A["Min.",1:3]),max(A["97.5%",1:3])),
		pch=15, cex=1.75, col="red",
		ylab="Population size", xlab="Estimation function", 
		bty="l", xaxt="n", las=1, font.lab=2, font=2, main="Means, medians & 95% quantile intervals")
	# Add median values to the plot
		points(x=c(1,2,3), y=A["Median",1:3], pch=16, cex=1.75, col="blue")
	# Draw the x-axis
		axis(1, c(1,2,3), tcl=0.3, labels=c("Kohn","Eggert","Chessel"), font=2)
	# Draw the 95% quantile interval
		arrows(c(1,2,3), A["2.5%",1:3], c(1,2,3), A["97.5%",1:3], length=0.11, angle=90, code=3, lwd=2)
	# Add legend
	legend("topright", c("Mean", "Median"), col=c(2,4), pch=c(15,16), cex=1, bty="n")
	}
