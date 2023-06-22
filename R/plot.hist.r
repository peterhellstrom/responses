#' @export
plot.hist <- function(obj, model, breaks=30, plot.type="histogram", drawN=TRUE) {

	rs <- obj$Resamplings
	nS <- obj$SampleSize
	ests <- obj$Coefficients[,1:3]
	colnames(ests) <- c("Kohn","Eggert","Chessel")

	# Plot distribution of estimates

	if (plot.type == "histogram") {
	# Histogram
		store.hist <- hist(ests[,model], freq=F, font.lab=2, las=1,
					breaks=breaks, col="lightgray", xlab="Estimate", ylab="Density",
					main=paste("Population size estimate,", model))
		if (drawN == TRUE) {
			# Overlay normal distribution
			.x <- seq(min(store.hist$breaks), max(store.hist$breaks), 0.1)
			lines(.x, dnorm(.x, mean=mean(ests[,model]), sd=sd(ests[,model])), type="l", lwd=2, col=4)
			}
		}
	if (plot.type == "density") {
		plot(density(ests[,model]), lwd=2, col=4, las=1, font.lab=2, main=paste("Population size estimate,", model))
		}

	# Vertical lines for max, min & median estimates
	abline(v=quantile(ests[,model],0.025), lty=2, lwd=2, col="red")
	abline(v=quantile(ests[,model],0.975), lty=2, lwd=2, col="red")
	abline(v=median(ests[,model]), lty=1, lwd=2, col="black")

	# Text showing mean value
	title(sub=paste("Median =", round(median(ests[,model]),2), ",", rs, "resamplings"), cex=0.8, font=3)
	intrvl <- round(as.numeric(quantile(ests[,model], c(0.025, 0.975))),2)
	mtext(text = paste("95% Quantile interval =", intrvl[1],"-",intrvl[2]), side=3, line=0.5)
}
