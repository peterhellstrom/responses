##################################################
# Install & load required packages
is.installed <- function(mypkg) mypkg %in% installed.packages()[,1]
install.missing <- function(mypkg, repos="http://ftp.sunet.se/pub/lang/CRAN/", dependencies=TRUE) {
	for (i in 1:length(mypkg)) if (is.installed(mypkg[i]) == FALSE) install.packages(pkgs=mypkg[i], lib=.Library, repos=repos, dependencies=dependencies)
}

pkgs <- c("MASS","fitdistrplus","fdrtool","R.utils","gamlss","gamlss.tr","bbmle","VGAM","xlsReadWrite")
install.missing(pkgs)

# Load required packages
library(MASS)
library(fitdistrplus)
library(fdrtool)
library(R.utils)
library(gamlss)
library(gamlss.tr)
library(bbmle)
library(VGAM)
library(xlsx)

##################################################
lnorm2norm <- function(mu,sigma) {
	names(mu) <- names(sigma) <- NULL
	 m <- exp(mu + (sigma^2)/2)
	 v <- exp(2*mu + sigma^2)*(exp(sigma^2)-1)
	 sd <- sqrt(v)
	 c(m=m, sd=sd, v=v)
}

norm2lnorm <- function(m,sd) {
	names(m) <- names(sd) <- NULL
	v <- sd^2
	mu <- log((m^2) / sqrt(v + m^2))
	sigma <- sqrt(log(v/(m^2) + 1))
	c(mu=mu,sigma=sigma,v=v)
}

##################################################
# Input is a list fitted with objects fitted with fitdist from fitdistrplus
# The input list must be created manually
# This function is an intermediate function, and calculates various information criteria
# It's intended use is ONLY from within extraDist, see below

extrDistrAIC <- function(x) {
	
	if (class(x)!="list") stop("Input must be a list")
	
	sapply(1:length(x), function(i) {
	
		n <- length(x[[i]]$data)
		k <- length(x[[i]]$estimate)

		loglik <- x[[i]]$loglik
	
		aic <- 2*k - 2*loglik
		# aic <- x[[i]]$aic
		aic.c <- -2*loglik + 2*k*(k+1)/(n-k-1)
		bic <- -2*loglik + k*log(n)
	
		.x <- c(n,k,loglik,aic,aic.c,bic)
	})
}
##################################################
# Calculates IC and IC weights, option order=TRUE/FALSE.
# order=TRUE returns the table with candidate models, in order according to AICc (from smallest to largest).

extrDistr <- function(x, order=FALSE) {
	.x <- extrDistrAIC(x)
	distname <- sapply(1:length(x), function(i) x[[i]]$distname)
	colnames(.x) <- distname
	rownames(.x) <- c("n","k","logLik","AIC","AICc","BIC")
	.x <- t(.x)
	#.x
	delta <- .x[,"AICc"] - min(.x[,"AICc"])
	L <- exp(-delta/2) # likelihoods of models
	w <-  L / sum(L) # Akaike weights
	.x <- data.frame(.x,delta,w)
	
	if (order==TRUE) {
		inds <- order(delta)
		.x <- .x[inds,]
		}
	.x
}
##################################################
# Plot cumulative density function
# Possible to plot several distributions in the same plot, by using a list:

ecdfn <- function(x) {

if (class(x) == "matrix") {
	sortx <- sapply(1:ncol(x), function(i) sort(x[,i]))
	ecdfx <- sapply(1:ncol(x), function(i) (1:length(x[,i])) / length(x[,i]))
	
	
	plot(sortx, ecdfx, type="n", xlab="x", ylab="Fn(x)")
	abline(h=c(0,1), lty=2, col="grey")
	if (!is.null(colnames(x))) legend("bottomright", colnames(x), pch=16, cex=1, col=1:length(x), bty="n")
	points(sortx, ecdfx, pch=16, cex=0.7, col=rep(1:ncol(x),each=nrow(x)))
	for (i in 1:ncol(x)) lines(sortx[,i], ecdfx[,i], col=i)
}

if (class(x) == "list") {
	sortx <- lapply(1:length(x), function(i) sort(x[[i]]))
	ecdfx <- lapply(1:length(x), function(i) (1:length(x[[i]])) / length(x[[i]]))
	
	
	plot(unlist(sortx), unlist(ecdfx), type="n", xlab="x", ylab="Fn(x)")
	abline(h=c(0,1), lty=2, col="grey")
	if (!is.null(names(x))) legend("bottomright", names(x), pch=16, cex=1, col=1:length(x), bty="n")
	for (i in 1:length(x)) points(sortx[[i]], ecdfx[[i]], pch=16, cex=0.7, col=i)
	for (i in 1:length(x)) lines(sortx[[i]], ecdfx[[i]], col=i)
	}

	out <- list(x=x, sortx=sortx, ecdfx=ecdfx)
	out
}
##################################################
# Another function for plotting an ecdf
ecdf2 <- function(x, plot=TRUE) {
	sortx <- sort(x)
	ecdfx <- (1:length(x)) / length(x)
	if (plot == TRUE) {
		
			plot(sortx, ecdfx, type="n", xlab="x", ylab="Fn(x)")
			abline(h=c(0,1), lty=2, col="grey")
			lines(sortx, ecdfx)
			points(sortx, ecdfx, pch=16, cex=.7)
		}
	out <- list(x=x,sortx=sortx,ecdfx=ecdfx)
	out
	}
##################################################
# Generate random variates from a truncated three-parameter t-distribution
rt3tr <- function(n, spec="t", a, b, mu, sigma, nu) {
	if (spec != "t") stop("You can only use the t-distribution")
	# Standardize truncation points
	a.st <- (a - mu) / sigma
	b.st <- (b - mu) / sigma
	# Generate random variates
	x <- (rtrunc(n=n, spec="t", df=nu, a=a.st, b=b.st) * sigma) + mu
	x
}

# PDF for three-parameter t-distribution (NOT truncated)
dt3 <- function(x,mu,sigma,nu) {
	lambda <- 1/sigma^2 # squared inverse scale parameter (precision)
	(gamma((nu+1)/2) / gamma(nu/2)) * ((lambda / (pi*nu))^(1/2)) * ((1 + ((lambda*(x-mu)^2) / nu))^(-(nu+1)/2))
}

##################################################
# Generate random variates from a theoretical distribution.
# distr = name of function, follow conventional R methods. E.g. "norm", "pois", "lnorm", "exp", etc.
# The function uses the random generation function for each distribution (if it exists), e.g. for distr="norm", "rnorm" is used.
# para = named list with parameters, must match arguments of the specified distribution.
# E.g. for distr="norm", see get("rnorm",mode="function") or args("rnorm")
# Allows sampling from groups, specified with the argument g.
# If a vector g is not supplied, the function uses g and prob to generate a random sample.
# Also possible with weighted sampling, argument prob (length of vector must equal number of categories in g).
# Output options: print=TRUE (prints data frame with grouping variable and random variates), print=FALSE (only random variates).

# Example:
#.x <- genRandom(distr="norm", para=list(n=10000, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, print=TRUE)
#hist(.x$x, col="lightgrey", breaks=60)

# A lot of the code is borrowed from the package fitdistrplus

genRandom <- function(distr, para, catgs=3, g=NULL, prob=NULL, print=FALSE, ...) {
	
	# Evaluate input arguments
	if(class(para) != "list") stop("'para' must be a named list")
	if (is.null(names(para))) stop("'para' must be a named list")
	if ((missing(distr) & !missing(para)) || (missing(distr) & !missing(para))) stop("'distr' and 'para' must defined")
	
	#if (!is.null(g)) if((class(g) %in% c("factor","integer")==FALSE)) stop("Grouping variable 'g' must be numeric")
	
	# Random generation from the distribution specified in 'distr'
	# & Check that the specified function exists.
	rdistname <- paste("r",distr,sep="")
	
	if (distr != "trunc") {
		if (!exists(rdistname, mode = "function")) stop(paste("The ", rdistname, " function must be defined"))
		# Check that the input parameters in the list 'para' matches the internal arguments of the function specified.
		# If not, stop and return error message.
		densfun <- get(rdistname, mode = "function")
		nm <- names(para)
		f <- formals(densfun)
		args <- names(f)
		m <- match(nm, args)
		if (any(is.na(m))) stop(paste("'para' specifies names which are not arguments to ", rdistname))
	}
	
	# If vector with probabilities is supplied, check that the number of categories matches length of prob. vector
	if (!missing(prob)) {
		if (length(prob) != catgs) stop("Length of 'prob' vector must equal number of categories")
	}
	
	# Create grouping variable, only if a grouping variable is not supplied by the user:
	if (is.null(g)) {
		if (missing(catgs)) stop("'catgs', number of categories must be entered")
		g <- sample(1:catgs, para[["n"]], replace=TRUE, prob=prob)
	}
	
	# For truncated distributions
	if (distr == "trunc" | distr == "t3tr") {	
		
		# Add code for case when spec is missing in rt3tr
		#if (distr == "t3tr") ... insert(para, 2, c(spec="t"))
		
		para.names <- names(para)
		para.temp <- lapply(3:length(para), function(i) para[[i]][g])
		para.temp <- insert(para.temp, ats=1, values=c(para$n), useNames=FALSE)
		para.temp <- insert(para.temp, ats=2, values=c(para$spec), useNames=FALSE)
		names(para.temp) <- para.names
		para <- para.temp
		
		x <- do.call(paste("r",distr,sep=""), para)
	}
	
	# For non-truncated distributions
	else if (distr != "trunc" | distr != "t3tr") {
		para.names <- names(para)
		para <- lapply(1:length(para), function(i) para[[i]][g])
		names(para) <- para.names
		x <- do.call(rdistname, para)
	}
	else stop("Not a valid distribution")
	
	# If print=TRUE, return a data frame with the grouping variable and not only the random variates
	if (print == TRUE) x <- data.frame(g,x)
	# Print
	x
}

##################################################
# Plot the density function of a theoretical distribution and also an area under the curve.
# distr = name of distribution, see documentation for function genRandom
# para = named list with parameters, must match arguments of distr
# xlim = vector of length 2, plot the density function from xlim[1] to xlim[2].
# from, to = plot area under curve in the region 'from' --> 'to'

# Example:
# auc("norm", para=list(mean=0, sd=1), from=qnorm(0.975, mean=0, sd=1), to=6, by=1/1000)

auc <- function(distr, para, xlim=c(-3,3), from, to, by=0.01, col="skyblue") {
	
	# Evaluate input arguments
	if(class(para) != "list") stop("'para' must be a named list")
	if (is.null(names(para))) stop("'para' must be a named list")
	if ((missing(distr) & !missing(para)) || (missing(distr) & !missing(para))) stop("'distr' and 'para' must defined")
	if (is.null(xlim)) stop("no valid x values, change the 'xlim' argument")
	
	# & Check that the specified function exists.
	ddistname <- paste("d",distr,sep="")
	if (!exists(ddistname, mode = "function")) stop(paste("The ", ddistname, " function must be defined"))
	
	# Check that the input parameters in the list 'para' matches the internal arguments of the function specified.
	# If not, stop and return error message.
	densfun <- get(ddistname, mode = "function")
	nm <- names(para)
	f <- formals(densfun)
	args <- names(f)
	m <- match(nm, args)
	if (any(is.na(m))) stop(paste("'para' specifies names which are not arguments to ", ddistname))
	
	# Create plot
	para.poly <- para
	x.poly <- seq(from,to,by)
	para.poly$x <- x.poly
	y.poly <- do.call(ddistname, para.poly)
	cord.x <- c(from, x.poly, to)
	cord.y <- c(0, y.poly, 0)

	x <- seq(xlim[1],xlim[2],by)
	para$x <- x
	y <- do.call(ddistname, para)
	plot(x,y,type="n",ylab="Density",las=1,main=distr)
	polygon(cord.x, cord.y, col=col, lty=2)
	lines(x,y,lty=1,col=1)
	
}

auc.poly <- function(distr,para,from,to,by=0.01) {
	# Evaluate input arguments
	if(class(para) != "list") stop("'para' must be a named list")
	if (is.null(names(para))) stop("'para' must be a named list")
	if ((missing(distr) & !missing(para)) || (missing(distr) & !missing(para))) stop("'distr' and 'para' must defined")
	
	# & Check that the specified function exists.
	ddistname <- paste("d",distr,sep="")
	if (!exists(ddistname, mode = "function")) stop(paste("The ", ddistname, " function must be defined"))
	
	# Check that the input parameters in the list 'para' matches the internal arguments of the function specified.
	# If not, stop and return error message.
	densfun <- get(ddistname, mode = "function")
	nm <- names(para)
	f <- formals(densfun)
	args <- names(f)
	m <- match(nm, args)
	if (any(is.na(m))) stop(paste("'para' specifies names which are not arguments to ", ddistname))
	
	para.poly <- para
	x.poly <- seq(from,to,by)
	para.poly$x <- x.poly
	y.poly <- do.call(ddistname, para.poly)
	cord.x <- c(from, x.poly, to)
	cord.y <- c(0, y.poly, 0)
	
	cbind(x=cord.x, y=cord.y)
	
}

##################################################
# CODE FROM Nadarajah & Kotz 2006 j stat soft

# this function computes (1) for given x, spec, a & b
dtrunc <- function(x, spec, a = -Inf, b = Inf, ...) {
    tt <- rep(0, length(x))
    g <- get(paste("d", spec, sep = ""), mode = "function")
    G <- get(paste("p", spec, sep = ""), mode = "function")
    tt[x>=a & x<=b] <- g(x[x>=a&x<=b], ...)/(G(b, ...) - G(a, ...))
    return(tt)
}

# this function computes (2) for given spec, a & b
extrunc <- function(spec, a = -Inf, b = Inf,...) {
    f <- function(x) x * dtrunc(x, spec, a = a, b = b, ...)
    return(integrate(f, lower = a, upper = b)$value)
}

# this function computes (3) for given spec, a & b
vartrunc <- function(spec, a = -Inf, b = Inf, ...) {
    ex <- extrunc(spec, a = a, b = b, ...)
    f <- function(x) (x - ex)^2 * dtrunc(x, spec, a = a, b = b, ...)
    tt <- integrate(f, lower = a, upper = b)$value
    return(tt)
}

# this function computes (4) for given x, spec, a & b
ptrunc <- function(x, spec, a=-Inf, b=Inf, ...) {
    tt <- x
    aa <- rep(a, length(x))
    bb <- rep(b, length(x))
    G <- get(paste("p", spec, sep = ""), mode = "function")
    tt <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
    tt <- tt - G(aa, ...)
    tt <- tt/(G(bb, ...) - G(aa, ...))
    return(tt)
}

# this function computes (5) for given spec, a & b
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
    tt <- p
    G <- get(paste("p", spec, sep = ""), mode = "function")
    Gin <- get(paste("q", spec, sep = ""), mode = "function")
    tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
    return(tt)
}

# this function computes (6) for given n, spec, a & b
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
    x <- u <- runif(n, min = 0, max = 1)
    x <- qtrunc(u, spec, a = a, b = b,...)
    return(x)
}
##################################################
##################################################

# Extract coefficients from a gamlss object
# Check help for link functions for the various parameters (necessary to transform to get the correct estimates)
# Only a few distributions supported, written for a very specific purpose.

get.coef <- function(object, method="coef") {
	pars.fitted <- object$parameters
	npar <- length(pars.fitted)
	if (method %in% c("coef","link") == FALSE) stop("Method must be 'coef' or 'link'")
	if (method == "coef") str <- paste("as.numeric(fitted(object,'", pars.fitted, "')[1])", sep="")
	if (method == "link") str <- rep(paste("object$", pars.fitted, ".link", sep=""))

	pars <- sapply(1:npar, function(i) { 
		val <- eval(parse(text=str[i]))
		ifelse(!is.null(val), val, NA)
		})
	names(pars) <- pars.fitted
	pars
}

coef.Gamlss <- function(object) {
	family <- object$family[1]
	pars.names <- object$parameters
	
	pars <- get.coef(object, "coef")
	pars.link <- get.coef(object, "link")
	npar <- length(pars.names); pars.tr <- numeric(npar)
	
	# Links: identity, log, logit. Not implemented yet: probit, cloglog
	invlogit <- function(x) exp(x)/(1 + exp(x))
	logit <- function(x) log(x / (1-x))
	
	if ( any(na.omit(pars.link) %in% c("identity","log","logit") == FALSE)) stop("Only 'identity', 'log' and 'logit' link functions are implemented yet")
	
	for (i in 1:npar) {
		pars.tr[i] <- switch(pars.link[i],
			'identity' = pars[i],
			'log' = log(pars[i]),
			'logit' = logit(pars[i]),
			'NA' = pars[i]) }
	
	names(pars) <- names(pars.tr) <- names(pars.link) <- pars.names
	data.frame(estimate=pars.tr, transformed=pars, link=pars.link)
	
}

##################################################
# Plot a gamlss.tr object and the data

plot.Gamlss <- function(object, a=NULL, b=NULL, breaks=50) {
	
	x <- as.numeric(object$y)
	
	if (is.null(a)) a <- min(x)
	if (is.null(b)) b <- max(x)
	
	n <- length(x)
	distr <- object$family[1]
	pars.names <- rownames(coef.Gamlss(object))
	pars <- as.numeric(coef.Gamlss(object)$transformed)
	names(pars) <- pars.names

	pars.str <- paste(pars.names,"=",pars)
	pars.str <- paste(rep(pars.str),collapse=", ")

	qstr <- paste("q",distr,"(p=ppoints(n),", pars.str, ")", sep="")
	qx <- eval(parse(text=qstr))
	qx.resid <- residuals(object) # ?residuals.gamlss, GAMLSS returns quantile residuals
	qx.resid <- qx.resid[!is.infinite(qx.resid)] # Remove infinite values

	dstr <- paste("d",distr,"(x,", pars.str, ")", sep="")
	dstr <- paste("curve(expr=", dstr, ", from=a, to=b, n=1001, add=T, col=2, lwd=2)" ,sep="")

	par(mfrow=c(2,2))
		hist(x, col="steelblue", breaks=breaks, freq=FALSE, main=paste("f:",distr))
		lines(density(x))
		eval(parse(text=dstr))
		
		plot(qx, sort(x), xlab="Theoretical quantiles", ylab="Sample quantiles", main=paste("qq:",distr), cex=0.75)
		abline(0,1,lty=2)
		
		hist(qx.resid, col="steelblue", breaks=breaks, freq=FALSE, xlab="Quantile residuals", main=paste("resid(",distr, ")", sep=""))
		curve(dnorm(x), n=1001, col=2, lwd=2, add=T)
		
		qqnorm(qx.resid, xlab="Theoretical quantiles", ylab="Sample quantiles", main=paste("qq: resid(",distr, ")", sep=""), cex=0.75)
		qqline(qx.resid, lty=2)
	par(mfrow=c(1,1))

}

##################################################
# Compare actual observed data with another (simulated) data set.
# Draw histogram of simulated data, and add min, max, density of observed data.
evalSim.plot <- function(x, x.sim) {
	
	x.eval.max <- 100 * (table(x.sim >= max(x)) / length(x.sim)) # % of simulated values larger than max(observed)
	x.eval.min <- 100 * (table(x.sim <= min(x)) / length(x.sim)) # % of simulated values smaller than min(observed)
	
	hist(x.sim, breaks=30, col="steelblue", xlab="x", ylab="Density", main="Simulation vs. data", freq=FALSE)
		lines(density(x.sim), lwd=2)
		lines(density(x),col=3, lwd=2, lty=1)
		abline(v=min(x), col=4, lwd=2, lty=2)
		abline(v=max(x), col=2, lwd=2, lty=2)
		legend("topleft", c("Max","Min","Density"), title="Observations", col=c(2,4,3), lwd=2, lty=c(2,2,1), bty="n", cex=0.75)
	title(sub=paste(x.eval.min[2], "&", as.numeric(x.eval.max[2]),"% of simulated random variates are outside observed range"))
}
##################################################
# Plot untruncated, left- & right truncated & right-truncated data:
# Inputs mu and sigma must be list with names: ut, lr, % r.
# a is lower boundary, b upper boundary.

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

	dev.new(width=12, height=6)
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

##################################################
# END
##################################################

