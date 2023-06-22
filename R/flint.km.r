# Revision of older functions (surv.prob.xxx) for Flint et al's Kaplan-Meier type estimator
# Staggered entry and censoring

#' @export
flint.km.nc <- function(x) {
	# No censoring and no missing data
	# Add warning messages!
	z <- ncol(x)
	M <- nrow(x)
	names.list <- list('Brood' = rownames(x), 'Time' = paste(colnames(x)[1:z-1], colnames(x)[2:z], sep="-"))

	# Point estimation
	S <- as.matrix(x[,-1] / x[,-z], ncol=z-1)
	dimnames(S) <- names.list
	w <- as.matrix(x[,-z], ncol=z-1)
	dimnames(w) <- names.list
	w <- (w*S) / S # weights
	S.t <- colSums(w*S, na.rm=TRUE) / colSums(w, na.rm=TRUE)
	S.f <- cumprod(S.t)
	#n.bar <- colMeans(x)[-z]
	#SE.t <- sqrt(colSums(x[,-z]^2 * t((t(S)-S.t)^2)) / (M * n.bar^2 * (M - 1)))
	#SE.t <- sqrt(colSums(x[,-z]^2 * sweep(S,2,S.t,"-")^2) / (M * n.bar^2 * (M - 1)))

	out <- t(data.frame(S.t, S.f))
	out
}

#' @export
flint.km.nc.boot <- function(x, nboot=2000, ci.type="perc") {
	tmp <- boot(data=x, statistic=function(x,i) flint.km.nc(x[i,]), R=nboot, stype="i")
	dims <- dim(tmp$t0)
	x.bar <- matrix(apply(tmp$t,2,mean), nrow=dims[1], ncol=dims[2], dimnames=dimnames(tmp$t0))
	sd.bar <- matrix(apply(tmp$t,2,sd), nrow=dims[1], ncol=dims[2], dimnames=dimnames(tmp$t0))
	# Confidence intervals
	ci <- sapply(1:ncol(tmp$t), function(i) boot.ci(tmp, type=ci.type, index=i)$percent)[4:5,]
	dimnames(ci) <- list(c("CI.L","CI.U"), paste(rep(dimnames(tmp$t0)[[1]],ncol(x)-1), rep(dimnames(tmp$t0)[[2]],each=ncol(x)-2)))

	ci.array <- array(dim=c(2,ncol(x)-1,2))
	ci.array[,,1] <- ci[,grep(colnames(ci),pattern="S.t")]
	ci.array[,,2] <- ci[,grep(colnames(ci),pattern="S.f")]

	dimnames(ci.array) <- list(c("CI.L","CI.U"), dimnames(tmp$t0)[[2]], dimnames(tmp$t0)[[1]])
	list('x.bar' = x.bar, 'sd.bar' = sd.bar, 'ci' = ci.array)
}

flint.km <- function(x, nboot=2000, ci.type="perc") {
	t1 <- flint.km.nc(x)
	t2 <- flint.km.nc.boot(x, nboot=nboot, ci.type=ci.type)
	out <- list('point'=t1, 'interval'=t2)
	class(out) <- "flint"
	out
}

#' @export
plot.flint <- function(object) {
	# Plotting of SE should be added (dotted or hashed lines)
	# object is an object fitted with surv.prob()
	if (length(colnames(object$x)) == 0) stop ("Please submit colnames to frame x corresponding to time")
	# Point estimates
	x.val <- as.numeric(colnames(object$x))
	y.val <- c(1, object$S.f)

	# Lines connecting point estimates
	lx.val <- rep(x.val,each=2)
	lx.val <- lx.val[-1]

	ly.val <- rep(y.val,each=2)
	ly.val <- ly.val[-length(ly.val)]

	plot(x=x.val, y=y.val, type="n", font.lab=2, las=1, xlim=c(0,max(x.val)), ylim=range(y.val), xlab="Time", ylab="Survival",
		main = "Modifed Kaplan-Meier estimate")
	points(x.val,y.val,pch=16,cex=1.4)
	lines(lx.val,ly.val)
}
