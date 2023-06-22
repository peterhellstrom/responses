# Analyze numerical response and includes non-linear fit for direct response
# Three steps:
# 1) plots input data
# 2) analyzes direct numerical response (isocline)
# 3) analyzes indirect numerical response (demographic or growth rate)

# Finally, fitted objects and model selection tables is printed

# SOME CAVEATS:
# For GAM analyses, the parameter k is set to -1 if n >= 10. Arbitrarily set to 4 if n < 10.

#' @export
nr <- function(x,y,plot=TRUE,method="direct") {

	# Create vector to be used in plots
	xv <- seq(0,max(x),0.1)

	# Get length of response variable and set GAM parameter k:
	n <- length(y)
	if (n < 10) {k <- 4}
	if (n >= 10) {k <- -1}

	# STEP 1
	nr.gam <- mgcv::gam(y ~ s(x, fx=FALSE, k=k, bs="cr"))
	# nr.micmen <- getInitial(y ~ SSmicmen(x, a, b),data=data.frame(x,y))
	a <- max(y)
	b <- mean(x)

	if (plot==TRUE) {
		# Plot simulated data (prey, predator, auto-correlations, cross-correlation)
		op <- c(
			par(mfrow=c(2,2)),
			par(mar=c(5,4,3,1)))

		plot(x=c(0,max(x)), y=c(0,max(y)),type="n",bty="l",
			xlab="Food density",ylab="Predator density",font.lab=2,las=1,ylim=c(0,max(y)),main="Isocline numeric response")
		lines(xv,predict(nr.gam,list(x=xv),type="response"),lwd=2,col=2)
		points(x,y,bty="l",font.lab=2,las=1)

		plot(x=c(0,max(x)), y=c(0,max(y)),type="n",bty="l",
			xlab="Food density",ylab="Predator density",font.lab=2,las=1,ylim=c(0,max(y)),main="Isocline numeric response")
		s <- 1:length(y)
		arrows(x[s],y[s],x[s+1],y[s+1], length=0.075, lty=1)

		plot(x,type="l",xlab="Time",ylab="Food density",bty="l",font.lab=2,las=1, main="Prey dynamics")
		plot(y,type="l",xlab="Time",ylab="Predator density",bty="l",font.lab=2,las=1, main="Predator dynamics")

		par(mfrow=c(1,1))
		par(op)

		op <- c(
			par(mfrow=c(3,2)),
			par(mar=c(5,4,4,1)))

		acf(x, main="ACF, prey", font.lab=2, las=1)
		pacf(x, main="PACF, prey", font.lab=2, las=1)
		acf(y, main="ACF, predator", font.lab=2, las=1)
		pacf(y, main="PACF, predator", font.lab=2, las=1)
		ccf(x,y,main="CCF, predator & prey", font.lab=2, las=1)

		par(mfrow=c(1,1))
		par(op)
	}

	# STEP 2
	# Direct numerical response
	nobs <- length(y)
	inf <- which(is.infinite(log(x)))
	if(length(inf) > 0) print("Warning: x variable contain zeroes, log-values infinite")

	fm0 <- lm(y ~ 1)
	fm1 <- lm(y ~ x)
	fm2 <- try(lm(y ~ log(x)))
	fm3 <- lm(y ~ x-1)
	fm4 <- try(nls(y ~ type2(x,a,b), start=list(a=a,b=b)))
	fm5 <- try(nls(y ~ type2t(x,a,b,c), start=list(a=a,b=b,c=min(x))))
	fm6 <- try(mgcv::gam(y ~ s(x, fx=FALSE, k=k, bs="cr")))

	out.dnr <- list(fm0,fm1,fm2,fm3,fm4,fm5,fm6)
	names(out.dnr) <- c("fm0","fm1","fm2","fm3","fm4","fm5","fm6")

	out.classes <- unlist(sapply(1:length(out.dnr), function(i) class((out.dnr)[[i]])))[1:7]
	out.inds <- which(out.classes!="try-error")

	test <- out.dnr[out.inds]
	modsel.dnr <- ICtab(test,delta=TRUE,weights=TRUE,sort=FALSE,type="AICc",nobs=nobs,mnames=names(test))

	if (plot==TRUE) {
		par(mfrow=c(2,4))

		plot(x,y, font.lab=2, las=1, xlab="Prey density", ylab="Predator density", main="fm0")
		abline(h=coef(fm0))

		plot(x,y, font.lab=2, las=1, xlab="Prey density", ylab="Predator density", main="fm1")
		abline(fm1)

		try(plot(log(x),y, font.lab=2, las=1, xlab="log(Prey density)", ylab="Predator density", main="fm2"))
		try(abline(fm2))

		plot(x,y, font.lab=2, las=1, xlab="Prey density", ylab="Predator density", main="fm3")
		try(abline(fm3))

		plot(x,y, font.lab=2, las=1, xlab="Prey density", ylab="Predator density", main="fm4")
		try(lines(xv,predict(fm4,list(x=xv))))

		plot(x,y, font.lab=2, las=1, xlab="Prey density", ylab="Predator density", main="fm5")
		try(lines(xv,predict(fm5,list(x=xv))))

		plot(x,y, font.lab=2, las=1, xlab="Prey density", ylab="Predator density", main="fm6")
		try(lines(xv,predict(fm6,list(x=xv),type="response")))

		par(mfrow=c(1,1))
	}

	if (method=="both") {
		# STEP 3
		# Indirect numerical response

		# Check if series contain zeros:
		x0 <- length(which(x==0))
		y0 <- length(which(y==0))

		if (x0 > 0) stop("x contains zeros, growth rates can't be calculated")
		if (y0 > 0) stop("y contains zeros, growth rates can't be calculated")

		# Calculate various variables
		# Growth rates
		ry <- diff(log(y))
		rx <- diff(log(x))

		ryc <- diff(y)
		rxc <- diff(x)

		# Non-lagged (t) and one-year lag (t1)
		yt <- y[-1]
		yt1 <- y[-length(y)]

		xt <- x[-1]
		xt1 <- x[-length(x)]

		# Ratio lagged prey/predator densities
		# ratxt1yt1 <- log(xt1)/log(yt1)
		ratxt1yt1 <- xt1/yt1
		ratio.warning <- length(which(is.infinite(ratxt1yt1)))
		if (ratio.warning > 0) stop("Ratio is infinite")

		# Residuals from a model with predator density vs. prey density
		resyt1 <- lm(y ~ x)
		resyt1 <- residuals(resyt1)[-length(fitted(resyt1))]

		logresyt1 <- lm(y ~ log(x))
		logresyt1 <- residuals(logresyt1)[-length(fitted(logresyt1))]

		nobs <- length(ry)

		fm1 <- lm(ry ~ xt)
		fm2 <- lm(ry ~ xt1)
		fm3 <- lm(ry ~ log(xt))
		fm4 <- lm(ry ~ log(xt1))
		fm5 <- lm(ry ~ ratxt1yt1)
		fm6 <- lm(ry ~ yt1)
		fm7 <- lm(ry ~ log(yt1))
		fm8 <<- lm(ry ~ xt + yt1)
		fm9 <- lm(ry ~ rx)
		fm9b <- lm(ryc ~ rxc)
		fm10 <<- lm(ry ~ rx + resyt1)
		fm11 <<- lm(ry ~ rx + logresyt1)

		out.inr <- list(fm1,fm2,fm3,fm4,fm5,fm6,fm7,fm8,fm9,fm9b,fm10,fm11)
		names(out.inr) <- c("fm1","fm2","fm3","fm4","fm5","fm6","fm7","fm8","fm9","fm10","fm11")

		modsel.inr <- ICtab(fm1,fm2,fm3,fm4,fm5,fm6,fm7,fm8,fm9,fm9b,fm10,fm11,delta=TRUE,weights=TRUE,sort=FALSE,type="AICc",nobs=nobs)
		r2.inr <- as.numeric(t(sapply(1:length(out.inr), function(i) summary(out.inr[[i]])$r.squared)))

		if (plot==TRUE) {
			op <- par(mfrow=c(2,5))

			plot(xt,ry, font.lab=2, las=1, xlab=expression(bold(paste("Prey density")[" t"])), ylab="Predator growth rate", main="fm1")
			abline(fm1)

			plot(xt1,ry, font.lab=2, las=1, xlab=expression(bold(paste("Prey density")[" t-1"])), ylab="Predator growth rate", main="fm2")
			abline(fm2)

			plot(log(xt),ry, font.lab=2, las=1, xlab=expression(bold(paste("log(Prey density)")[" t"])), ylab="Predator growth rate", main="fm3")
			abline(fm3)

			plot(log(xt1),ry, font.lab=2, las=1, xlab=expression(bold(paste("log(Prey density)")[" t-1"])), ylab="Predator growth rate", main="fm4")
			abline(fm4)

			plot(ratxt1yt1,ry, font.lab=2, las=1, xlab=expression(bold(paste("Prey density")[" t-1"] / paste("Predator density")[" t-1"])), ylab="Predator growth rate", main="fm5")
			abline(fm5)

			plot(yt1,ry, font.lab=2, las=1, xlab=expression(bold(paste("Predator density")[" t-1"])), ylab="Predator growth rate", main="fm6")
			abline(fm6)

			plot(log(yt1),ry, font.lab=2, las=1, xlab=expression(bold(paste("log(Predator density)")[" t-1"])), ylab="Predator growth rate", main="fm7")
			abline(fm7)

			plot(rx,ry, font.lab=2, las=1, xlab=expression(bold(paste("Prey growth rate"))), ylab="Predator growth rate", main="fm9")
			abline(fm9)

			plot(rxc,ryc, font.lab=2, las=1, xlab=expression(bold(paste("Prey change"))), ylab="Predator change", main="fm9b")
			abline(fm9b)

			par(mfrow=c(1,1))

			par(op)

			# Plot response with two response variables, 3d-planes
			par(mfrow=c(1,3))

			curve3d(
				expr=coef(fm8)[1] + coef(fm8)[2]*x + coef(fm8)[3]*y,
				from = c(min(xt),min(yt1)), to = c(max(xt),max(yt1)), n = c(25,25),
				theta = -50, phi = 30,
				xlab = "Prey density t", ylab = "Predator density t-1", zlab="Predator growth rate",
				font.lab = 2, las = 1,
				main = "fm8, Prey- and predator-dependent indirect nr",
				sys3d = "persp",
				box = TRUE,
				axes = TRUE,
				ticktype = "detailed"
			)

			curve3d(
				expr=coef(fm10)[1] + coef(fm10)[2]*x + coef(fm10)[3]*y,
				from = c(min(rx),min(resyt1)), to = c(max(rx),max(resyt1)), n = c(25,25),
				theta = -50, phi = 30,
				xlab = "Prey growth rate", ylab = "Residual predator density", zlab="Predator growth rate",
				font.lab = 2, las = 1,
				main = "fm10, Prey- and predator-dependent indirect nr",
				sys3d = "persp",
				box = TRUE,
				axes = TRUE,
				ticktype = "detailed"
			)

			curve3d(
				expr=coef(fm11)[1] + coef(fm11)[2]*x + coef(fm11)[3]*y,
				from = c(min(rx),min(logresyt1)), to = c(max(rx),max(logresyt1)), n = c(25,25),
				theta = -50, phi = 30,
				xlab = "Prey growth rate", ylab = "Log(Residual predator density)", zlab="Predator growth rate",
				font.lab = 2, las = 1,
				main = "fm11, Prey- and predator-dependent indirect nr",
				sys3d = "persp",
				box = TRUE,
				axes = TRUE,
				ticktype = "detailed"
			)

		par(mfrow=c(1,1))
		}
	}

	# Output for direct-model only
	if (method=="direct") {
		return(list(
			# Fitted models
			fm.dnr = out.dnr,
			# Model selection
			aicc.dnr = modsel.dnr
		))}

	if(method=="both") {
		# Output
		return(list(
			# Fitted models
			fm.dnr = out.dnr,
			fm.inr = out.inr,
			# Model selection
			aicc.dnr = modsel.dnr,
			aicc.inr = modsel.inr,
			r2.inr = data.frame(model = names(out.inr), r2 = round(r2.inr,3))
		))}

}
