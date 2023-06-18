# Plot correlation between parameter estimates for Kohn and Eggert models:
plot.par.cor <- function(obj, ...) {
	dev.new(width=12,height=6)
	par(mfrow=c(1,2))
	plot(obj$Coefficients[,"b.Kohn"], obj$Coefficients[,"a.Kohn"], cex=0.7, 
	xlab="Half-saturation constant", ylab="Asymptote", font.lab=2, las=1, main="Kohn")
	abline(h=obj$ParameterSummary["Mean","a.Kohn"], col=2, lty=2, ...)

	plot(obj$Coefficients[,"b.Eggert"], obj$Coefficients[,"a.Eggert"], cex=0.7,
	xlab="Natural logarithm of the rate constant", ylab="Asymptote", font.lab=2, las=1, main="Eggert")
	abline(h=obj$ParameterSummary["Mean","a.Eggert"], col=2, lty=2, ...)
	par(mfrow=c(1,1))
}
