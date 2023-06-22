# ANALYZE DENSITY-DEPENDENCE IN PREDATION RATE
# Fit polynomial regression of order 0-3 to the data.

# General function
#' @export
pred.rate <- function(x, y, plot=TRUE) {

	x <- x
	p <- y/x

	fit0 <- lm(p ~ 1)
	#fit1 <- lm(p ~ x)
	#fit2 <- lm(p ~ x + I(x^2))
	#fit3 <- lm(p ~ x + I(x^2) + I(x^3))

	fit1 <- lm(p ~ poly(x, degree=1, raw=TRUE))
	fit2 <- lm(p ~ poly(x, degree=2, raw=TRUE))
	fit3 <- lm(p ~ poly(x, degree=3, raw=TRUE))

	print(list(
		fit0 = summary(fit0),
		fit1 = summary(fit1),
		fit2 = summary(fit2),
		fit3 = summary(fit3)
	))

	print(anova(fit0,fit1,fit2,fit3))

	if (plot) {
		xp <- seq(0,max(x),0.1)
		yp1 <- predict(fit1, data.frame(x=xp))

		yp2 <- predict(fit2, data.frame(x=xp))

		yp3 <- predict(fit3, data.frame(x=xp))

		plot(x,p,ylab="Relative predation rate",pch=16,font.lab=2,las=1,
			main="Polynomial regression: predation rate",
			ylim=range(as.vector(c(yp1,yp2,yp3))))
		lines(xp,yp1,lty=1)
		lines(xp,yp2,lty=2)
		lines(xp,yp3,lty=3)
		legend("topright", c("1st order","2nd order","3rd order"),lty=c(1,2,3),bty="n")
	}
}

#degree <- 3
#strFn <- "x + "
#for (i in 2:degree) strFn <- paste(strFn, "I(x^",i,") + ",sep="")
