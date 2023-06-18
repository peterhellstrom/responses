# The two functions below provide plots of simple analyses of covariance,
# with only one factor and one covariate. Both functions plot color coded
# data points from the different groups, and the first also plots fitted
# with separate slopes, whereas for the second there is a common slope

plot.ancova <- function(response, covariate, g, model=c("1","2","3"), cex=1.5, lwd=2, dev.new=TRUE) {

	model <- match.arg(model)
	
	if (dev.new) dev.new()
	plot(covariate, response, cex=cex, lwd=lwd, col=as.numeric(g)+1)
	grps <- unique(g)
	ng <- length(grps)
	if (model == "1") {
		# Separate slopes and intercept
		for (i in 1:ng) abline(lm(response[g==grps[i]] ~ covariate[g==grps[i]]), lwd=2, col=as.numeric(grps[i])+1)
	} else if (model == "2") {
		# Separate slopes without intercept
		for (i in 1:ng) abline(lm(response[g==grps[i]] ~ covariate[g==grps[i]] - 1), lwd=2, col=as.numeric(grps[i])+1)
	} else if (model == "3") {
		# Common slopes, different intercept
	}
}

# common slope plot
plot.ancova.com <- function(response,covariate,g) {
	g <- as.factor(g)
	dev.new()
	plot(covariate, response, cex=1.5, lwd=2, col=as.numeric(g)+1)
	grps <- unique(g)
	ng <- length(grps)
	fit.dat <- data.frame(resp=response, grp=g, cov=covariate)
	fm <- lm(resp ~ grp + cov, data=fit.dat)
	slp <- coef(fm)[ng+1]
	pred.dat <- data.frame(grp=levels(g), cov=rep(0.0, ng))
	intr <- predict(fm, pred.dat)
	for (i in 1:ng)
		abline(intr[i], slp, lwd=2, col=as.numeric(grps[i])+1)
}
