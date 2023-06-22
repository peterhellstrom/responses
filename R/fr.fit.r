# Fit functional responses:
# Input arguments:
# data: list containing list objects called x & y
# method: nls or mle
# eq: parameterization, Michaelis-Menten or Holling
# Output: x, y, coefficients, and fitted objects are stored in output
# NOTE: Model selection was previously performed within fr.fit, but was moved
# to an external function, aicc.c, that can be used to get a table with AICc, model ranks etc.
# mle estimation is performed with the BFGS

# + See fr.fit.selfStart.r

# UPDATED 2013-11-14
# Major revision of maximum likelihood functions
# sigma must now be included in the estimation procedure, start values are required also for sigma
# see corresponding updates in fr.fit.selfStart.r
#' @export
fr.fit <- function(x, y, method=c("nls","mle"), eq=c("M-M","Holling"), nls.start,
	iter=5000, tol=0.000001, minFactor=0.000001) {

xy <- data.frame(x=x,y=y)
method <- match.arg(method)
eq <- match.arg(eq)

# Control iterations etc in nls:
nlc <- nls.control(maxiter = iter, tol=tol, minFactor=minFactor)

if (eq == "M-M") str.colnames <- c("a","b","theta","sigma")
else if (eq == "Holling") str.colnames <- c("a","h","theta","sigma")

	if (method=="nls") {
		if (eq=="M-M") {
			fm.0 <- lm(y~1)
			fm.I <- lm(y~x-1)
			fm.II <- nls(y ~ SStypeII.nls(x,a,b), data=xy, control=nlc)
			fm.IIIa <- nls(y ~ SStypeIIIa.nls(x,a,b), data=xy, control=nlc)
			fm.IIIb <- nls(y ~ SStypeIIIb.nls(x,a,b,theta), data=xy, control=nlc)
		} else if (eq=="Holling") {
			fm.0 <- lm(y~1)
			fm.I <- lm(y~x-1)
			fm.II <- nls(y ~ TypeII.h(x,a,h), start=nls.start[["ini.II"]][1:2], data=xy, control=nlc)
			fm.IIIa <- nls(y ~ TypeIIIa.h(x,a,h), start=nls.start[["ini.IIIa"]][1:2], data=xy, control=nlc)
			fm.IIIb <- nls(y ~ TypeIIIb.h(x,a,h,theta), start=nls.start[["ini.IIIb"]][1:3], data=xy, control=nlc)
		}

		fit.objects <- list(fm.0=fm.0, fm.I=fm.I, fm.II=fm.II, fm.IIIa=fm.IIIa, fm.IIIb=fm.IIIb)
		p.est <- sapply(fit.objects, coef)
		names(p.est[[1]]) <- "a"
		names(p.est[[2]]) <- "a"
		coefs <- with(p.est, rbind.fill.matrix(t(fm.0), t(fm.I), t(fm.II),t(fm.IIIa), t(fm.IIIb)))
		coefs[3:4,3] <- 1:2

	} else if (method=="mle") {
		if (eq=="M-M") {
			# Obtain mle's
			fm.0 <- mle2(minuslogl=LL.null, start=as.list(nls.start[["ini.0"]]), data=as.list(xy), control=list(maxit=iter))
			fm.I <- mle2(minuslogl=LL.l, start=as.list(nls.start[["ini.I"]]), data=as.list(xy), control=list(maxit=iter))
			fm.II <- mle2(minuslogl=LL.mm, start=as.list(nls.start[["ini.II"]]), fixed=list(theta=1), data=as.list(xy), control=list(maxit=iter))
			fm.IIIa <- mle2(minuslogl=LL.mm, start=as.list(nls.start[["ini.IIIa"]]), fixed=list(theta=2), data=as.list(xy), control=list(maxit=iter))
			fm.IIIb <- mle2(minuslogl=LL.mm, start=as.list(nls.start[["ini.IIIb"]]), data=as.list(xy), control=list(maxit=iter))
		} else if (eq=="Holling") {
			fm.0 <- mle2(minuslogl=LL.null, start=as.list(nls.start[["ini.0"]]), data=as.list(xy), control=list(maxit=iter))
			fm.I <- mle2(minuslogl=LL.l, start=as.list(nls.start[["ini.I"]]), data=as.list(xy), control=list(maxit=iter))
			fm.II <- mle2(minuslogl=LL.h, start=as.list(nls.start[["ini.II"]]), fixed=list(theta=1), data=as.list(xy), control=list(maxit=iter))
			fm.IIIa <- mle2(minuslogl=LL.h, start=as.list(nls.start[["ini.IIIa"]]), fixed=list(theta=2), data=as.list(xy), control=list(maxit=iter))
			fm.IIIb <- mle2(minuslogl=LL.h, start=as.list(nls.start[["ini.IIIb"]]), data=as.list(xy), control=list(maxit=iter))
		}

		fit.objects <- list(fm.0=fm.0, fm.I=fm.I, fm.II=fm.II, fm.IIIa=fm.IIIa, fm.IIIb=fm.IIIb)
		p.est <- sapply(fit.objects, coef)
		coefs <- with(p.est, rbind.fill.matrix(t(fm.0), t(fm.I), t(fm.II),t(fm.IIIa), t(fm.IIIb)))
		coefs <- coefs[,c(1,3,4,2)]
	}

	# Output
	rownames(coefs) <- c("0","I","II","IIIa","IIIb")

	out <- list(method=method, eq=eq, x=x, y=y, coefs=coefs, models=fit.objects)
	class(out) <- "fResponse"
	return(out)
}
