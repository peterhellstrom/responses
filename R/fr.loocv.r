# Cross validation & jack-knifing
# Leave one out cross validation = loocv

# x = prey density
# y = kill rate
# starts = start estimates for nls, list
# type = type of functional response, Holling type II or Holling type IIIa or IIIb

fr.loocv <- function(x, y, names=NULL, starts, type=c("TypeII","TypeIIIa","TypeIIIb"), plot=TRUE, plot.parameters=c("a","b")) {

	type <- match.arg(type)
	# Store output
	err <- matrix(ncol=9, nrow=length(y))
	colnames(err) <- c("Deviance", "x.obs", "y.obs", "y.pred", "PredictionError", "y.pred.full", "a", "b", "theta")
	rownames(err) <- names(x)

	if (type=="TypeII") {
		total.fit <- nls(y ~ TypeII(x,a,b), start=list(a=starts[1],b=starts[2]))
		parms.total <- c(coef(total.fit), "TotalDeviance"=deviance(total.fit))
	}
	if (type=="TypeIIIa"){
		total.fit <- nls(y ~ TypeIIIa(x,a,b), start=list(a=starts[1],b=starts[2]))
		parms.total <- c(coef(total.fit), "TotalDeviance"=deviance(total.fit))
	}
	if (type=="TypeIIIb"){
		total.fit <- nls(y ~ TypeIIIb(x,a,b,theta), start=list(a=starts[1],b=starts[2],theta=starts[3]))
		parms.total <- c(coef(total.fit), "TotalDeviance"=deviance(total.fit))
	}

	for (i in 1:length(y)) {

		# Remove each data point for each loop step, "Training" data set:
		y.train <- y[-i]
		x.train <- x[-i]

		# Fit the desired functional response with nls
		if (type=="TypeII") {
			fit <- nls(y.train ~ TypeII(x.train,a,b), start=list(a=starts[1],b=starts[2]))
			coefs <- c(coef(fit),NA)
		} else if (type=="TypeIIIa") { 
			fit <- nls(y.train ~ TypeIIIa(x.train,a,b), start=list(a=starts[1],b=starts[2]))
			coefs <- c(coef(fit),NA)
		} else if (type=="TypeIIIb") { 
			fit <- nls(y.train ~ TypeIIIb(x.train,a,b,theta), start=list(a=starts[1],b=starts[2],theta=starts[3]))
			coefs <- coef(fit)
		}    

		# Calculate output
		# 1) Residual sum of squares when data point i is removed
		err[i,1] <- deviance(fit) # equals sum(resid(fit)^2)
		# 2) Squared difference between observed and predicted y
		# The removed data point's x value (prey density)
		x.obs <- x[i]
		# Take the difference between the observed yi at xi, subtract the predicted yi
		# (from fit) and square this difference
		y.obs <- y[i]
		y.pred <- predict(fit, data.frame(x.train = x.obs))
		y.pred.full <- fitted(total.fit)[i]
		err[i,2] <- x.obs
		err[i,3] <- y.obs
		err[i,4] <- y.pred 
		err[i,5] <- (y.obs - y.pred)^2
		err[i,6] <- y.pred.full
		err[i,7:9] <- coefs
	}

	# Print output
	out <- list(
		x = x, y=y,
		fr.type = type, loocv = data.frame(err), mean.deviance = mean(err[,1]),
		mean.prediction.error = mean(err[,2]), sd.deviance = sd(err[,1]),
		sd.prediction.error = sd(err[,2]), total = parms.total
	)
	
	if (is.null(names)) names <- 1:length(y)
	
	if (plot) {
		dev.new(width=8,height=5)
		par(mfrow=c(1,2))
		
		plot(x, y, 
			xlim=c(0, max(x)), ylim=c(0,max(c(y,err[,4]))), 
			font.lab=2, las=1, pch=16, font.lab=2, las=1, main="Predictions")
		points(err[,"x.obs"], err[,"y.pred"], pch=16,col=2)
		legend("bottomright", c("Observed", "Jack. predictions"), col=c(1,2), pch=c(16, 16), bty="n", cex=0.8)
	
		plot(err[,plot.parameters[1]], err[,plot.parameters[2]], pch=16, font.lab=2,
			xlab=plot.parameters[1], ylab=plot.parameters[2], main="Jack-knife influence plot")
		text(err[,"a"], err[,"b"], labels = names, pos=1, offset=0.5, cex=0.7)
		points(parms.total[1], parms.total[2], pch=16, cex=1.6, col=2)
	
		par(mfrow=c(1,1))
	}
	
	out
}
