# Function that calculates confidence interval around a (linear) regression line
# and outputs vectors that can be plotted as shaded polygons
	
line.ci <- function(model, xlims, by=0.1) {
	xv <- seq(xlims[1],xlims[2],by)
	variable <- attr(model$terms,"term.labels")
	xv <- list(xv)
	names(xv) <- variable
	predci <- predict(model,xv,interval=c("confidence"),level=0.95)
	CI.U <- predci[, "upr"] 
	CI.L <- predci[, "lwr"] 
	xv <- data.frame(xv)
	names(xv) <- "x"
	xvec <- c(xv$x, tail(xv$x, 1), rev(xv$x), xv$x[1]) 
	yvec <- c(CI.L, tail(CI.U, 1), rev(CI.U), CI.L[1])
	out <- list(x=xvec,y=yvec,predci=predci)
	out
}

# Same as line.ci, but normal approximation to non-linear function
line.ci.nls <- function(model, xlims, by=0.1) {
	if (class(model) != "nls") stop("Only nls class allowed")

	xv <- seq(from=xlims[1], to=xlims[2], by=by)
	variable <- attr(model$dataClasses, "names")
	xv <- list(xv)
	names(xv) <- variable
	se.fit <- sqrt(apply(attr(predict(model,xv),"gradient"), 1, function(x) sum(vcov(model) * outer(x,x))))

	predci <- predict(model, xv) + outer(se.fit,qnorm(c(.5, .025,.975)))
	colnames(predci) <- c("fit","lwr","upr")
	CI.U <- predci[,2] 
	CI.L <- predci[,3]
	xv <- data.frame(xv)
	names(xv) <- "x"
	xvec <- c(xv$x, tail(xv$x, 1), rev(xv$x), xv$x[1]) 
	yvec <- c(CI.L, tail(CI.U, 1), rev(CI.U), CI.L[1])
	out <- list(x=xvec, y=yvec, predci=predci)
	out
}

lm.plot <- function(x) {
	dev.new()
	par(mfrow=c(2,2))
	plot(x)
	par(mfrow=c(1,1))
}
