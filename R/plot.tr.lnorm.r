# Plot untruncated, left- & right truncated & right-truncated data:
# Inputs mu and sigma must be a list with names: ut, lr, % r.
# a is lower boundary, b upper boundary.
#' @export
plot.tr.lnorm <- function(mu, sigma, a, b) {

	if (class(mu)!="list") stop("Input must be a list")
	if (is.null(names(mu))) stop("'mu' must be a named list")
	if (length(mu) != 3) stop("Input list must be of length=3")

	if (class(sigma)!="list") stop("Input must be a list")
	if (is.null(names(sigma))) stop("'sigma' must be a named list")
	if (length(sigma) != 3) stop("Input list must be of length=3")

	xv <- seq(0,b,0.01)
	yv.nt <- dnorm(xv, mean=mu$ut, sd=sigma$ut)
	yv.rl <- dtrunc(xv, spec="lnorm", a=a, b=b, meanlog=mu$lr, sdlog=sigma$lr)
	yv.r <- dtrunc(xv, spec="lnorm", a=-Inf, b=b, meanlog=mu$r, sdlog=sigma$r)
	ylims <- c(0, max(yv.nt,yv.rl,yv.r))

	par(mfrow=c(1,2))
	curve(dlnorm(x, meanlog=mu$ut, sdlog=sigma$ut), n=1001, col=1, add=F, from=0, to=b, ylim=ylims, ylab="Density", main="Log-normal distribution")
	curve(dtrunc(x, "lnorm", a=a, b=b, meanlog=mu$lr, sdlog=sigma$lr), n=1001, col=2, lty=2, add=T)
	curve(dtrunc(x, "lnorm", a=-Inf, b=b, meanlog=mu$r, sdlog=sigma$r), n=1001, col=3, lty=3, add=T)
	legend("topright", c("Not truncated", "Left- & right truncated", "Right truncated"), col=1:3, lty=1:3, bty="n", title="Truncation",cex=0.8)

	# on log-scale

	xlim.log <- c(qnorm(0.001,mu$ut,sigma$ut),qnorm(0.999,mu$ut,sigma$ut))
	xv.log <- seq(xlim.log[1],xlim.log[2],0.01)
	yv.nt.log <- dnorm(xv.log, mean=mu$ut, sd=sigma$ut)
	yv.rl.log <- dtrunc(xv.log, spec="norm", a=log(a), b=log(b), mean=mu$lr, sd=sigma$lr)
	yv.r.log <- dtrunc(xv.log, spec="norm", a=-Inf, b=log(b), mean=mu$r, sd=sigma$r)
	ylims.log <- c(0, max(yv.nt.log,yv.rl.log,yv.r.log))

	curve(dnorm(x, mean=mu$ut, sd=sigma$ut), n=1001, col=1, add=F, xlim=xlim.log, ylim=ylims.log, ylab="Density", xlab="log(x)", main="Log-normal distribution, log-scale")
	curve(dtrunc(x, "norm", a=log(a), b=log(b), mean=mu$lr, sd=sigma$lr), n=1001, col=2, lty=2, add=T)
	curve(dtrunc(x, "norm", a=-Inf, b=log(b), mean=mu$r, sd=sigma$r), n=1001, col=3, lty=3, add=T)
	legend("topright", c("Not truncated", "Left- & right truncated", "Right truncated"), col=1:3, lty=1:3, bty="n", title="Truncation",cex=0.8)
	par(mfrow=c(1,1))
}
