
########################################################################
# Functions for numerical response analyses
# Updated 2014-03-31



########################################################################
# Required packages for use of nResponse
# Check that required packages are installed

#is.installed <- function(mypkg) mypkg %in% installed.packages()[,1]
#install.missing <- function(mypkg, repos="http://ftp.acc.umu.se/mirror/CRAN/", dependencies=TRUE) {
#	for (i in 1:length(mypkg)) if (is.installed(mypkg[i]) == FALSE) install.packages(pkgs=mypkg[i], lib=.Library, repos=repos, dependencies=dependencies)
#}

# Exclude glmmadmb from list of packages (it's not available from CRAN For R 3.x)
#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")), type="source")

#pkgs <- c("bbmle","emdbook","mgcv","boot","simpleboot","lme4","lmerTest","MCMCglmm","nnet","xlsx","nlme","coefplot2","pscl","RLRsim","gstat","sp","plotrix","lattice","RODBC","ltm")
#is.installed(pkgs)
#install.missing(pkgs)

# Load required packages
#pkgs <- as.list(c(pkgs, "glmmADMB"))
#lapply(pkgs, require, character.only=TRUE)



########################################################################
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



########################################################################
pn <- function(x) {crossprod(!is.na(x))}

cor.prob <- function(x, method=c("pearson")) {
	method <- match.arg(method)
	# Correlations below main diagonal
	# Significance tests with pairwise deletion above main diagonal
	pair.SampSize <- pn(x)
	above1 <- row(pair.SampSize) < col(pair.SampSize)
	pair.df <- pair.SampSize[above1] - 2
	R <- cor(x, use="pair", method)
	above2 <- row(R) < col(R)
	if (method=="pearson") {
		r2 <- R[above2]^2
		Fstat <- (r2 * pair.df)/(1 - r2)
		R[above2] <- 1 - pf(Fstat, 1, pair.df)
	}
	R
}

# correlation function 
# based on post by Bill Venables on R-Help 
# Date: Tue, 04 Jan 2000 15:05:39 +1000 
# https://stat.ethz.ch/pipermail/r-help/2000-January/009758.html
# modified by G L Simpson, September 2003 
# version 0.2: added print.cor.prob 
#              added class statement to cor.prob 
# version 0.1: original function of Bill Venables 
corProb <- function(X, dfr = nrow(X) - 2) { 
    R <- cor(X) 
    above <- row(R) < col(R) 
    r2 <- R[above]^2 
    Fstat <- r2 * dfr / (1 - r2) 
    R[above] <- 1 - pf(Fstat, 1, dfr) 
    class(R) <- "corProb" 
    R 
}

print.corProb <- function(x, digits = getOption("digits"), quote = FALSE, na.print = "", 
    justify = "none", ...) { 
    xx <- format(unclass(round(x, digits = 4)), digits = digits, justify = justify) 
    if (any(ina <- is.na(x))) 
        xx[ina] <- na.print 
    cat("\nCorrelations are shown below the diagonal\n") 
    cat("P-values are shown above the diagonal\n\n") 
    print(xx, quote = quote, ...) 
    invisible(x) 
} 



########################################################################
# Combined standardized indices

# Citations:
# de la Mare, W.K., and Constable, A.J. 2000. Utilising data from ecosystem monitoring for managing fisheries: 
# development of statistical summaries of indices arising from the CCAMLR ecosystem monitoring program. CCAMLR Sci. 7: 101-117.
# Boyd, I.L., and Murray, A.W.A. 2001. Monitoring a marine ecosystem using responses of upper trophic level predators. J. Anim. Ecol. 70: 747-760.
# Reid, K., Croxall, J.P., Briggs, D.R., and Murphy, E.J. 2005. 
# Antarctic ecosystem monitoring: quantifying the response of ecosystem indicators to variability in Antarctic krill. ICES J. Mar. Sci. 62: 366-373.
# See also chapter 11 by Croxall in
# Boyd, I., Wanless, S., and Camphuysen, C.J. 2006. Top Predators in Marine Ecosystems.

csi <- function(x) {
	# From Box 11.3 in Top Predators in Marine Ecosystems, eds. I.L.Boyd,S. Wanless and C. J. Camphuysen.
	# Croxall: Monitoring predator-prey interactions using multiple predator species: the South Georgia experience
	# Missing values? How to deal with them?
	
	x <- scale(x)
	a <- ifelse(!is.na(x), 1, 0)
	
	I <- sapply(1:nrow(x), function(i) t(a[i,]) %*% x[i,])
	# Edit: check that matrix is positive and semi-definite (some eigenvalues are negative)
	S <- cov(x)
	V <- diag(a %*% S %*% t(a))
	csi <- I / sqrt(V)
	list(csi=csi, I=I, V=V, S=S, eig=eigen(S), x=x)
}



########################################################################
# Inverse logit
expit <- function(x) exp(x) / (1+exp(x))


########################################################################
# Revision of older functions (surv.prob.xxx) for Flint et al's Kaplan-Meier type estimator
# Staggered entry and censoring

flint.km.nc <- function(x) {
	# No censoring and no missing data
	# Add warning messages!
	z <- ncol(x)
	M <- nrow(x)
	names.list <- list('Brood' = rownames(x), 'Time' = paste(colnames(x)[1:z-1], colnames(x)[2:z], sep="-"))
	
	# Point estimation
	S <- as.matrix(x[,-1] / x[,-z], ncol=z-1); dimnames(S) <- names.list
	w <- as.matrix(x[,-z], ncol=z-1); dimnames(w) <- names.list
	w <- (w*S) / S # weights
	S.t <- colSums(w*S, na.rm=TRUE) / colSums(w, na.rm=TRUE)
	S.f <- cumprod(S.t)
	#n.bar <- colMeans(x)[-z]
	#SE.t <- sqrt(colSums(x[,-z]^2 * t((t(S)-S.t)^2)) / (M * n.bar^2 * (M - 1)))
	#SE.t <- sqrt(colSums(x[,-z]^2 * sweep(S,2,S.t,"-")^2) / (M * n.bar^2 * (M - 1)))
	
	out <- t(data.frame(S.t, S.f))
	out
}

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



########################################################################
# NOTE! This function does not account for irregular census intervals
# E.g. time series must have only a single deltat / frequency!
gr <- function(N, method="r", plot=TRUE) {

	if (class(N) != "ts") stop("Object N must be a time series object")
	
	na.inds <- which(is.na(N))
	if (length(na.inds) > 0) message("Object N contains missing/NA values")
	
	g <- switch(method,
		lambda = N[-1]/N[-length(N)],
		r = diff(log(N))
	)
	
	# Time component of the variable g is wrong (it is associated with t instead of t-1)
	# and must therefore be changed
	g <- ts(g, start=start(N)[1], deltat=deltat(N))
	
	title.text <- switch(method,
		lambda = expression(bold(paste("Growth rate =", lambda))),
		r = expression(bold(paste("Growth rate = log(", lambda, ")")))
	)
		
	if (plot) {
	
		# LOESS for N
		loess.N <- loess(N ~ time(N))
		# LOESS for g
		loess.g <- loess(g ~ time(g))
		# LOESS for dd in g
		dd.N <- N[-length(N)]
		loess.dd <- loess(g ~ dd.N)
		pred.dd.N <- seq(min(N,na.rm=T),max(N,na.rm=T),length=101)
		pred.loess.dd <- predict(loess.dd, newdata=pred.dd.N)
		
		# Density dependence in growth rates
		lm.dd <- lm(g ~ N[-length(N)])
		coefs <- summary(lm.dd)$coefficients[,1]
		p <- summary(lm.dd)$coefficients[2,4]
		lm.text <- paste(round(coefs[1],6), " + ", round(coefs[2],6), "x, p =", round(p,6), sep="")
		
		dev.new()
		par(mfrow=c(2,2))
		plot(N, las=1, type="n", main="Population counts", xlab="Time", ylab="Counts")
			lines(N, lty=2)
			points(N, pch=16, cex=1)
			points(time(N), predict(loess.N), type="l", col=2, lwd=2)
		plot(g,las=1,type="n",main=title.text,xlab="Time",ylab="Growth rate")
			lines(g,lty=2)
			points(g,pch=16,cex=1)
			points(time(g), predict(loess.g), type="l", col=2, lwd=2)
		plot(N[-length(N)],las=1,g,pch=16,cex=1,xlab="Population density",ylab="Growth rate",main="Density-dependence in growth rates")
			points(pred.dd.N,pred.loess.dd,col=2,lwd=2,type="l")
		acf(g, main="Auto-correlation in growth rates", na.action=na.pass)
		if (length(na.inds) > 0) title(sub="N contains NA values")
		par(mfrow=c(1,1))
	}
	# Print output
	invisible(list(TimeSeries=N, GrowthRate=g))
}



########################################################################
list2ascii <- function(x,file=paste(deparse(substitute(x)),".txt",sep="")) { 

   tmp.wid = getOption("width")  # save current width
   options(width=10000)          # increase output width
   sink(file)                    # redirect output to file
   print(x)                      # print the object
   sink()                        # cancel redirection
   options(width=tmp.wid)        # restore linewidth
   return(invisible(NULL))       # return (nothing) from function
}

listToArray <- function(L, dimnames=NULL) {
	a <- array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
	dimnames(a) <- dimnames
	a
}

tapply.formula <- function(fo, df, func=mean, output=c("matrix", "data.frame"), rpl.NA=FALSE) {
	# fo = formula
	# df = data.frame
	# func = function
	
	output <- match.arg(output)
	
    mf <- model.frame(fo, df)
    i <- attr(attr(mf, 'terms'), 'response')
    y <- mf[[i]]
    y.name <- colnames(mf)[i]
    by <- mf[-i]

    # return(as.data.frame.table(tapply(y, by, func, na.rm=TRUE), responseName=y.name))
	
	if (output == "data.frame") {
		out <- as.data.frame.table(tapply(y, by, func), responseName=y.name)
	} else if (output == "matrix") {
		out <- tapply(y, by, func)
	}
	
	if (rpl.NA != FALSE) out[is.na(out)] <- 0
	out
}


# http://stackoverflow.com/questions/7196450/create-a-data-frame-of-unequal-lengths

cfun <- function(L) {
	pad.na <- function(x, len) {
		c(x,rep(NA,len-length(x)))
	}
	maxlen <- max(sapply(L, length))
	do.call(data.frame,lapply(L, pad.na, len=maxlen))
}

na.pad <- function(x,len){
    x[1:len]
}

makePaddedDataFrame <- function(l,...){
    maxlen <- max(sapply(l, length))
    data.frame(lapply(l, na.pad, len=maxlen),...)
}

# Example:
# L <- list(x = c(rep("one", 2)), y = c(rep("two", 10)), z = c(rep("three", 5)))
# makePaddedDataFrame(L)
# t(cfun(L))

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

reorder2 <- function(x, X){
	X <- rev(X)
	for(i in seq_along(X)) x <- relevel(x, X[i])
	x
}

ICtab.df <- function(x) {
	if (class(x) != "ICtab") stop("Input must be an ICtab")
	z <- as.data.frame(matrix(unlist(x), ncol=length(x), nrow=length(x[[1]])))
	colnames(z) <- names(x)
	rownames(z) <- attr(x, "row.names")
	z
}


########################################################################
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


########################################################################
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




########################################################################
lm.plot <- function(x) {
	dev.new()
	par(mfrow=c(2,2))
	plot(x)
	par(mfrow=c(1,1))
}



########################################################################
# logit
logIt <- function(x) 1/(1+1/(exp(x)))



########################################################################
# Mayfield - estimator
# Flint et al 1995
#######################################################################

mayfield <- function(dat, method) {

	n.broods <- dim(dat)[1]
	cols <- ncol(dat)

	if (is.null(ncol(dat)) == TRUE) stop("Data frame must have at least two columns!")
	n.obs <- apply(dat, 1, function(x) length(which(!is.na(x))))
	n.intervals <- n.obs - 1

	# Calculate delta brood survival
	# Find all observations on a row
	inds <- lapply(1:n.broods, function(i) which(!is.na(dat[i,])))
	dat.obs <- lapply(1:n.broods, function(i) dat[i,inds[[i]]])

	# Calculate brood size at first observation - brood size at final observation
	delta.bs <- sapply(1:n.broods, function(i) as.numeric(dat.obs[[i]][1] - dat.obs[[i]][length(dat.obs[[i]])]))

	time <- unlist(strsplit(colnames(dat),"X"))
	time <- as.numeric(time[time!=""])
	time.obs <- lapply(1:n.broods, function(i) time[inds[[i]]])

	# If observations are made at equally spaced intervals:
	L <- sapply(1:length(time.obs), function(i) unique(diff(time.obs[[i]])))
	n.L <- length(unique(n.obs))

	###
	if (method == "exposure-days") {

		if (n.L == 1) {
			cat(paste("Broods were observed at equally spaced intervals, L = ", unique(L)), "\n")
			if (cols == 2) sr <- dat[,-1] / dat[,-ncol(dat)]
			if (cols > 2) sr <- rowSums(dat[,-1]) / rowSums(dat[,-ncol(dat)])
			dsr <- sr^(1/L)
			exposure <- delta.bs / (1-dsr)
			exposure2 <- rowSums(L*(dat[-1] + dat[-ncol(dat)])*0.5)
			# exposure not defined if delta.bs == 0, and dsr == 1. Returns NaN
			# Replace with values from exposure2
			nan.inds <- which(is.nan(exposure))
			exposure[nan.inds] <- exposure2[nan.inds]
		}

		if (n.L > 1) {
			cat("Broods were not observed at equally spaced intervals", "\n")
			# Calculate exposure for each brood
			# Step 1: Calculate sum of two consecutive counts / broods
			step1 <- lapply(1:n.broods, function(i) dat.obs[[i]][-length(dat.obs[[i]])] + dat.obs[[i]][-1])
			# Step 2: Calculate length of each interval
			step2 <- lapply(1:n.broods, function(i) diff(time.obs[[i]]))
			# Step 3: Calculate exposure (assuming changes in brood size occur at the midpoint of the observation interval)
			step3 <- lapply(1:n.broods, function(i) step1[[i]] * step2[[i]] * 0.5)
			# Final step, sum all exposure days
			exposure <- sapply(1:n.broods, function(i) sum(step3[[i]]))
			# Daily Survival Rate
			dsr <- 1 - (delta.bs / exposure)
		}

	} else if (method == "bart-robson") {

		pm.mayfield <- function(alive,dead,intrvl) 1 - (sum(dead) / sum((intrvl*(alive + dead*0.5))))
		dsr <- numeric(length(dat.obs))
		exposure <- numeric(length(dat.obs))

		for (i in 1:length(dat.obs)) {
			temp <- as.numeric(dat.obs[[i]])
			dead <- temp[-length(temp)] - temp[-1]
			alive <- temp[-length(temp)] - dead
			intrvl <- diff(time.obs[[i]])

			# Original estimate
			p.hat <- pm.mayfield(alive,dead,intrvl)
			if (p.hat == 1 & sum(dead)==0) {
				dsr[i] <- p.hat
			} else {
				n.iter <- 10
				p.iter <- numeric(n.iter)

				for (j in 1:n.iter) {
					f.pm <- sum((intrvl/p.hat) * (alive - (dead*p.hat^intrvl / (1 - p.hat^intrvl))))
					f.prim.pm <- sum((intrvl/p.hat^2) * (alive + ((dead*p.hat^intrvl)*(intrvl-1+p.hat^intrvl)) / (1-p.hat^intrvl)^2))
					p.hat <- p.hat + f.pm/f.prim.pm
					p.iter[j] <- p.hat
				} # end of inner loop
	
				dsr[i] <- p.iter[n.iter]
				# Problem here: exposure is not defined if dsr == 1, and prints as zero.
				exposure[i] <- delta.bs[i] / (1 - dsr[i])
			} # end of else

		} # end of outer loop

		# Calculate exposure2 for each brood
		# Step 1: Calculate sum of two consecutive counts / broods
		step1 <- lapply(1:n.broods, function(i) dat.obs[[i]][-length(dat.obs[[i]])] + dat.obs[[i]][-1])
		# Step 2: Calculate length of each interval
		step2 <- lapply(1:n.broods, function(i) diff(time.obs[[i]]))
		# Step 3: Calculate exposure (assuming changes in brood size occur at the midpoint of the observation interval)
		step3 <- lapply(1:n.broods, function(i) step1[[i]] * step2[[i]] * 0.5)
		# Final step, sum all exposure days
		exposure2 <- sapply(1:n.broods, function(i) sum(step3[[i]]))

		nan.inds <- which(exposure==0)
		exposure[nan.inds] <- exposure2[nan.inds]

	} # end of method bart-robson

	###
	# Calculate standard error
	se.dsr <- function(exposure, dsr) {
		M <- n.broods
		exposure.bar <- sum(exposure/M)
		dsr.bar <- 1 - (sum(delta.bs) / sum(exposure))
		num <- sum((exposure^2)*((dsr - dsr.bar)^2))
		denom <- ((exposure.bar^2)*(M - 1)*M)
		sqrt(num / denom)
	}

	dsr.pt <- 1 - (sum(delta.bs) / sum(exposure))
	dsr.int <- se.dsr(exposure,dsr)

	###
	# Collect output

	out <- list(
		# List with brood level output
		Brood = data.frame(n.obs, delta.bs, exposure, dsr),
		# Vector with population level output
		PointEstimate = c(n = sum(n.obs), delta.bs = sum(delta.bs), exposure = sum(exposure), dsr = dsr.pt),
		# Vector with standard error and confidence interval
		IntervalEstimate = c(se = dsr.int, cil = dsr.pt - 2*dsr.int, ciu = dsr.pt + 2*dsr.int)
	) # End of output list

	# Print output
	out

}



########################################################################
# "Standard" type II & type III models
type2 <- deriv(~ (a*x) / (b + x), c("a","b"), function(x,a,b) {})
type3 <- deriv(~ a*x^2 / (b^2 + x^2), c("a","b"), function(x,a,b) {})
type3b <- deriv(~ a*x^theta / (b^theta + x^theta), c("a","b","theta"), function(x,a,b,theta) {})

# Threshold models
type2ta <- deriv(~ (a*(x-d)) / ((b-d) + (x-d)), c("a","b","d"), function(x,a,b,d) {})
type2tb <- deriv(~ (a*(x-d)) / (b + (x-d)), c("a","b","d"), function(x,a,b,d) {})

type3ta <- deriv(~ (a*(x-d)^theta) / ((b-d)^theta + (x-d)^theta), c("a","b","d","theta"), function (x,a,b,d,theta) {} )
type3tb <- deriv(~ (a*(x-d)^theta) / ((b)^theta + (x-d)^theta), c("a","b","d","theta"), function (x,a,b,d,theta) {} )

type3tb <- function (x,a,b,d,theta) (a*(x-d)^theta) / ((b)^theta + (x-d)^theta)

# Different parameterizations of the logistic function used by Henden et al.
logist <- deriv(~ a/(1 + exp(-(b + d*x))), c("a","b","d"), function(x,a,b,d) {})
logist2 <- deriv(~ a/(1 + exp((b + d*-x))), c("a","b","d"), function(x,a,b,d) {})
logist3 <- deriv(~ a/(1 + exp((b - d*x))), c("a","b","d"), function(x,a,b,d) {})



########################################################################
nr.lm.g <- function(x, y, g) {

	if (length(levels(g)) > 2) stop("Current version only allows two levels in factor g")
	
	# Plot data
	if (min(x) < 0) {xlims <- c(min(x),max(x))}
	if (min(x) >= 0) {xlims <- c(0,max(x))}
	ylims <- c(0,max(y))
	gnum <- as.numeric(g)
	gnum[which(gnum==2)] <- 16

	# Fit the different models
	fm1 <- lm(y ~ x * g) # Varying intercepts and slopes
	fm2 <- lm(y ~ g:x) # Varying slopes, common intercept
	fm3 <- lm(y ~ x + g) # Varying intercepts, common slope
	fm4 <- lm(y ~ g:x - 1) # No intercept, varying slopes
	fm5 <- lm(y ~ x) # No effect of grouping variable, only of covariate

	dev.new(width=13, height=6)
	op <- par(list(mfrow=c(2,3), mar=c(5,4,3,1)))

	plot(x, y, pch=gnum, main="fm1: Varying intercepts and slopes", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=coef(fm1)[1], b=coef(fm1)[2], lty=1)
	abline(a=coef(fm1)[1] + coef(fm1)[3], b=coef(fm1)[2] + coef(fm1)[4], lty=2)

	plot(x, y, pch=gnum, main="fm2: Varying slopes, common intercept", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=coef(fm2)[1], b=coef(fm2)[2], lty=1)
	abline(a=coef(fm2)[1], b=coef(fm2)[3], lty=2)

	plot(x, y, pch=gnum, main="fm3: Varying intercepts, common slope", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=coef(fm3)[1], b=coef(fm3)[2], lty=1)
	abline(a=coef(fm3)[1]+coef(fm3)[3], b=coef(fm3)[2], lty=2)

	plot(x, y, pch=gnum, main="fm4: No intercept, varying slopes", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=0, b=coef(fm4)[1], lty=1)
	abline(a=0, b=coef(fm4)[2], lty=2)

	plot(x, y, pch=gnum, main="fm5: No effect of covariate", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(fm5)

	par(mfrow=c(1,1))
	par(op)

	out.nr.lm <- list(fm1,fm2,fm3,fm4,fm5)
	names(out.nr.lm) <- c("fm1","fm2","fm3","fm4","fm5")
	modsel.nr.lm <- ICtab(fm1,fm2,fm3,fm4,fm5,delta=TRUE,sort=TRUE,weights=TRUE,type="AICc",nobs=length(fitted(fm1)))
	r2.nr.lm <- as.numeric(t(sapply(1:length(out.nr.lm), function(i) summary(out.nr.lm[[i]])$r.squared)))

	ind <- attr(modsel.nr.lm,"row.names")[1]
	ind2 <- which(names(out.nr.lm)==ind)

	# Output
	list(
		fm.nr.lm = out.nr.lm,
		aicc.nr.lm = modsel.nr.lm,
		r2.nr.lm = r2.nr.lm,
		best = anova(out.nr.lm[[ind2]],fm5,test="F")
	)
}



########################################################################
# Analyze numerical response and includes non-linear fit for direct response
# Three steps:
# 1) plots input data
# 2) analyzes direct numerical response (isocline)
# 3) analyzes indirect numerical response (demographic or growth rate)

# Finally, fitted objects and model selection tables is printed

# SOME CAVEATS:
# For GAM analyses, the parameter k is set to -1 if n >= 10. Arbitrarily set to 4 if n < 10.

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
		dev.new(width=8, height=6)

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

		dev.new(width=6,height=8)

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
		dev.new(width=10,height=6)
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
			dev.new(width=22,height=12)
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
			dev.new(width=22,height=12)
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



########################################################################
# Test for overdispersion (works with lme4 and glmmadmb)
# http://glmm.wikidot.com/faq

overdisp <- function(model) {
	# number of variance parameters in an n-by-n variance-covariance matrix
	vpars <- function(m) {
		nrow(m)*(nrow(m)+1)/2
	}
	model.df <- sum(sapply(VarCorr(model),vpars)) + length(fixef(model))
	rdf <- nrow(model.frame(model)) - model.df
	if (class(model) != "glmmadmb") {
		rp <- residuals(model, type="pearson")
	} else {
		rp <- model$residuals
	}
	Pearson.chisq <- sum(rp^2)
	prat <- Pearson.chisq / rdf
	pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
	c(chisq=Pearson.chisq, ratio=prat, rdf=rdf, p=pval)
}

# From Zuur et al. 2013, p. 138 [A Beginner's Guide to GLM and GLMM with R]
overdisp.zuur <- function(model, data) {
	if (class(model) != "glmmadmb") {
		E1 <- resid(model, type="pearson")
	} else {
		E1 <- model$residuals
	}
	N <- nrow(data)
	p <- length(fixef(model)) + 1
	od <- sum(E1^2) / (N - p)
	od
}

# overdisp.zuur(model=m1.admb, data=dat)


########################################################################
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

	dev.new(width=8, height=6)
	plot(x=x.val, y=y.val, type="n", font.lab=2, las=1, xlim=c(0,max(x.val)), ylim=range(y.val), xlab="Time", ylab="Survival",
		main = "Modifed Kaplan-Meier estimate")
	points(x.val,y.val,pch=16,cex=1.4)
	lines(lx.val,ly.val)
}



########################################################################
# Custom functions for standard error & coefficient of variation
se <- function(x) sd(x) / sqrt(length(x))
cv <- function(x) sd(x) / mean(x)
ci.l <- function(x) {
	mean(x) - qnorm(0.975)*se(x)
}
ci.u <- function(x) {
	mean(x) + qnorm(0.975)*se(x)
}

# ToDO: Handling of NA's in tapply.formula?
pop.samp <- function (fo, df, rownames=TRUE, ci=FALSE) {
	
	df <- droplevels(df)
	z <- tapply.formula(fo, df, length)
	dims <- dim(z)
	out.names <- dimnames(z)
	
	if (ci == TRUE) {
		funcs <- c("sum", "length", "mean", "median", "sd", "se", "cv", "min", "max", "ci.l", "ci.u")
		out <- sapply(1:length(funcs), function (i) do.call(tapply.formula, list(fo=fo, df=df, func=get(funcs[i]))))
		colnames(out) <- c("Sum", "n", "Mean", "Median", "SD", "SE", "CV", "min", "max", "CIL", "CIU")
	} else if (ci == FALSE) {
		funcs <- c("sum", "length", "mean", "median", "sd", "se", "cv", "min", "max")
		out <- sapply(1:length(funcs), function (i) do.call(tapply.formula, list(fo=fo, df=df, func=get(funcs[i]))))
		colnames(out) <- c("Sum", "n", "Mean", "Median", "SD", "SE", "CV", "min", "max")
	}

	if (length(dims) > 1) {
		g <- expand.grid(out.names)
		g.row <- apply(g[,1:ncol(g)], 1, paste, collapse = ".")
		out <- data.frame(out, g)
		if (rownames) rownames(out) <- g.row
		out
	} else {
		out <- as.data.frame(out)
		out <- data.frame(out, out.names)
		if (rownames == FALSE) rownames(out) <- NULL
		out
	}
}



########################################################################
ProjectSurvival <- function(N0, T, st) {
	# N0 = initial number of individuals
	# T = total number of transitions ("years"-1)
	# st = time specific survival (vector of same length as T)
	if (length(st) != T) stop("st must be a vector with length T")
	out = matrix(ncol=T+1, nrow=length(N0))
	out[,1] <- N0
	for (t in 1:T) {
		N.surv = rbinom(n=N0, size=N0, prob=st[t])
		out[,t+1] <- N.surv
		N0 = N.surv
	}
	colnames(out) <- paste("t",1:ncol(out),sep="")
	rownames(out) <- paste("b",1:nrow(out),sep="")
	out
}



########################################################################
# Convert to repeated measures input:
# Input is a data frame or matrix where rows are broods and time are columns
rep.bn <- function(x) {
	if (length(colnames(x)) == 0) stop("Please supply column names for input matrix")
	brood <- factor(1:nrow(x))
	rownames(x) <- brood
	z <- ncol(x)
	temp <- lapply(2:z, function(i) data.frame(Alive=x[,i], Dead=x[,i-1]-x[,i], brood=brood, time=colnames(x)[i]))
	temp <- do.call("rbind", temp)
	out <- list(
		surv = cbind(Alive=temp$Alive, Dead=temp$Dead),
		brood = factor(temp$brood),
		timestep = factor(temp$time)
	)
	out
}

#inp <- rep.bn(x)


########################################################################
repr.histories <- function(x) {
	# Input x is a n-dimensional array created by table()
	x <- as.data.frame(x)
	x <- x[x$Freq > 0,]
	rownames(x) <- NULL
	x
}



########################################################################
# Plot reproduction data (yearly averages), based on data summarized with the pop.samp function
repr.plot <- function(x, y, g, error=NULL, n, xlims=range(x, na.rm=T), ylims=NULL,
	xlab="Rodent density index", ylab="Population average", main="",
	mar = c(5,5.5,1,2), cex=1.3, line=3.5,
	points=TRUE, pch = as.numeric(g), dev.new=FALSE, legend=TRUE, ...) {
	
	if (is.null(ylims)) {
		if (is.null(error)) {
			ylims <- c(0,ceiling(max(y, na.rm=T)))
		} else {
			ylims <- c(0,ceiling(max(y + error, na.rm=T)))
		}
	}
	
	if (dev.new == TRUE) dev.new(width=7, height=6)
	op <- par(mar=mar)
	plot(x, y, type="n", xlab=xlab, ylab="", main=main, font.lab=2, las=1, cex.axis=cex, cex.lab=cex, ylim=ylims, xlim=xlims, bty="l")
	mtext(text=ylab, side=2, line=line, cex=cex, font=2)
	
	if (!is.null(error)) segments(x, y - error, x, y + error)
	if (points) points(x, y, pch=pch, cex=cex, ...)
	if (legend) repr.plot.legend()
	par(op)
}


repr.plot.points <- function(x, y, pch, error=NULL, cex=1.3, ...) {
	if (!is.null(error)) segments(x, y - error, x, y + error)
	points(x, y, pch=pch, cex=cex, ...)
}

repr.plot.legend <- function(pos="bottomright", cex=1.3) legend(pos, c("1970s","2000s"), pch=c(16, 21), pt.bg=c("black","white"), cex=cex, lty=c(1,2), bty="n")



########################################################################
# Observed data - Calculate survival
# Input columns with number of individuals alive at each time step and grouping variable g
# g can only contain two levels!
surv.prob.2rand <- function(x, g, R=5000, plot=TRUE) {
	if (length(levels(g)) > 2) stop("It's only possible to compare two survival functions")

	# Test statistic (Pollock et al 1989)
	surv.0 <- surv.prob.n(x=x, g=g)
	S.t.0 <- t(sapply(surv.0, "[[", "S.t"))
	SE.t.0 <- t(sapply(surv.0, "[[", "SE.t"))
	# Calculate test-statistic
	# z-score for each individual conditional survival time
	z0.t <- diff(apply(S.t.0, 2, rev)) / colSums(SE.t.0)
	# z.t <- diff(S.t.0) / colSums(SE.t.0)
	d2.obs <- sum(z0.t^2)

	# Randomization procedure
	d2.rand <- numeric(R) # Store output

	for (i in 1:R) {
		# Resample the grouping vector:
		g.rand <- sample(g, replace=FALSE, size=length(g))
		# Calculate survival for randomized sample
		surv.rand <- surv.prob.n(x=x, g=g.rand)
		S.t.rand <- t(sapply(surv.rand, "[[", "S.t"))
		SE.t.rand <- t(sapply(surv.rand, "[[", "SE.t"))
		# Calculate test-statistic
		# z-score for each individual conditional survival time
		z.t <- diff(apply(S.t.rand, 2, rev)) / colSums(SE.t.rand)
		# z.t <- diff(S.t.rand) / colSums(SE.t.rand)
		d2 <- sum(z.t^2)
		d2.rand[i] <- d2
	}

	# Get p-value
	p.rand <- length(which(d2.rand >= d2.obs)) / R

	if (plot) {
		# Histogram of test statistic z
		hist(d2.rand)
		abline(v=d2.obs, lty=2, col=2, lwd=2)
	}

	#out <- list(d2.obs=d2.obs, d2.rand=d2.rand, p=p.rand)
	out <- list(d2.obs=d2.obs, p=p.rand)
	out
}



########################################################################
# Bootstrap
q.fun <- function(x, na.rm=TRUE) quantile(x, c(0.025,0.5,0.975), na.rm=na.rm) # Function for extracting quantiles
randomSample <- function(df, n=nrow(df), replace=TRUE) df[sample(nrow(df), n, replace=replace),]

surv.prob.boot <- function(x, nboot=2000, na.rm=TRUE) {
	
	# Resample data, generate matrix with index variables (columns are replicates, rows are broods)
	inds <- sapply(1:nboot, function(i) sample(1:nrow(x), replace=TRUE, size=nrow(x)))
	# Point estimation for resampled data
	est <- lapply(1:nboot, function(i) surv.prob(x = x[inds[,i],]))
	# Extract data
	# Extract survival estimates (S.t) and survival function (S.f)
	S.t <- t(sapply(est, "[[", "S.t"))
	S.f <- t(sapply(est, "[[", "S.f"))
	
	# Calculate quantiles
	if (dim(x)[2] > 2) {
		out.boot <- list(
			S.t = t(apply(S.t, 2, q.fun, na.rm=na.rm)), S.f = t(apply(S.f, 2, q.fun, na.rm=na.rm)),
			SE.S.t = apply(S.t, 2, sd, na.rm=na.rm), SE.S.f = apply(S.f, 2, sd, na.rm=na.rm),
			nboot = nboot)
	} else {
		out.boot <- list(
			S.t = t(apply(S.t, 1, q.fun, na.rm=na.rm)), S.f = t(apply(S.f, 1, q.fun, na.rm=na.rm)),
			SE.S.t = apply(S.t, 1, sd, na.rm=na.rm), SE.S.f = apply(S.f, 1, sd, na.rm=na.rm),
			nboot = nboot)
	}
	out.boot
}

surv.prob.boot.n <- function(x, g, nboot=2000, na.rm=TRUE) {
	x.g <- split(x, g)
	lapply(x.g, surv.prob.boot, nboot=nboot, na.rm=na.rm)
}



########################################################################
# surv.prob
# KAPLAN-MEIER product-type estimator
# See also flint.km.r for revised versions of some of the functions
# After Flint et al. 1995
# TO DO:
# Add more error handling

# Kaplan-Meier type: Calculations for standard errors: Should censored broods be included or not???

# Known errors in Kaplan-Meier type function:
# 1) Code breaks if all broods in a bootstrap sample are censored (should be uncommon, but happens for small samples)
# Some errors were caused by me trying to name columns in bootstrap output.
# Now bootstrap output are given without column names corresponding to time.
# When bootstrapping with few observations per time unit a common warning is: In qt(p, df, lower.tail, log.p) : NaNs produced

# GENERAL POINTS OF IMPORTANCE:
# - Non-independence among siblings
# - Hatching order (possibly) influences survival probability
# - Non-constant probability of survival over time, i.e. if one chick dies, the probability of survival for the remaining chicks may increase or decrease.
# - Success depends on number of "trials"?
# - Consider a death as a "removal" - will survival probabilities follow a hypergeometric or multinomial distribution?
# - beta-binomial distribution for survival probability?, see Bolker's book!

surv.prob <- function(x) {

	# x is a matrix with broods as rows and time as columns
	# z = number of time steps, S = survival (proportion), w = weights, numbers at risk
	# censored, non.censored
	# S.t = weighted mean of the estimates for individual broods, S.f = Kaplan-Meier survival function estimates of the proportion surviving from day 0, 
	# SE.it, SE.t, 
	# n.bar = Mean brood size at time t, M = number of broods (or marked females with broods)

	z <- ncol(x) # Number of time steps = z-1
	counts <- nrow(x) # Number of rows/broods
	
	names.list <- list('Brood' = rownames(x), 'Time' = paste(colnames(x)[1:z-1], colnames(x)[2:z], sep="-"))
	
	# Calculate survival 
	S <- sapply(1:(z-1), function(i) x[,i+1] / x[,i]); dimnames(S) <- names.list
	# Numbers at risk, used as weights.
	w <- cbind(x[,1:(z-1)]); dimnames(w) <- names.list
	
	# Logical, look for undefined survival estimates in each brood
	for (i in 2:z) undefined <- length(which(is.nan(S)))
		
		if (undefined > 0) {
			# Remove undefined values, using the any(), is.nan(), and drop functions
			x <- as.matrix(x[apply(S*w, 1, function(x)!any(is.nan(x))),, drop=F])
			w <- as.matrix(w[apply(S*w, 1, function(x)!any(is.nan(x))),, drop=F])
			S <- as.matrix(S[apply(S, 1, function(x)!any(is.nan(x))),, drop=F])
		}

		# CHANGE HERE!
		# If a brood only has a single observation (NA), survival can't be calculated and that brood should be removed:
		# nas <- apply(S, 1, function(x) length(which(is.na(x) & !is.nan(x)))) 
		nas <- apply(S, 1, function(x) length(which(is.na(x)))) # Count number of NA's on each row
		nas.inds <- which(nas == ncol(S))

		if (length(nas.inds > 0)) {
			x <- as.matrix(x[-nas.inds,])
			S <- as.matrix(S[-nas.inds,])
			w <- as.matrix(w[-nas.inds,])
		}

		valid.counts <- dim(S)[1] - sapply(1:ncol(S), function(i) length(which(is.na(S)[,i])))
		names(valid.counts) <- names.list$Time
		censored <- which(is.na(S*w))
		c.count <- length(censored)

		if (c.count == 0) {
			S.t <- if (z == 2) (sum(S*w, na.rm=T) / sum(w)) else (colSums(S*w, na.rm=T) / colSums(w, na.rm=T))
			S.f <- cumprod(S.t)
			SE.it <- w^2 * t((t(S)-S.t)^2)
			n.bar <- colMeans(as.matrix(x[,-c(z)]), na.rm=T)
			M <- valid.counts
			SE.denom <- M * n.bar[-z]^2 * (M-1)
			SE.t <- sqrt(colSums(SE.it, na.rm=T) / SE.denom[1:(z-1)])
			non.censored <- valid.counts
		} else if (c.count > 0) {
			if (z == 2) {
				S.t <- sum(S*w, na.rm=T) / sum(w[-censored])
				S.f <- cumprod(S.t)
				SE.it <- w^2 * t((t(S)-S.t)^2)
				# CODE STOPS HERE if bootstrap sample ends up with all samples being censored.
				n.bar <- colMeans(x[-censored,-z], na.rm=T)
				M <- valid.counts - length(censored)
				SE.denom <- M * n.bar[-z]^2 * (M-1)
				SE.t <- sqrt(colSums(SE.it, na.rm=T) / SE.denom[1:(z-1)])
				non.censored <- valid.counts - length(censored)
			} else if (z > 2) { 		
				w <- (S*w)/S
				S.t <- colSums(S*w, na.rm=T) / colSums(w, na.rm=T)
				S.f <- cumprod(S.t)
				SE.it <- w^2 * t((t(S)-S.t)^2)
				xS <- ifelse(!is.na(S),1,NA) # A matrix with zero's and ones indicating valid, non-censored observations (==1)
				M <- apply(xS,2,sum,na.rm=T)
				n.bar <- colMeans(xS * x[,-z], na.rm=T)
				SE.denom <- M * n.bar^2 * (M-1)
				SE.t <- sqrt(colSums(SE.it, na.rm=T) / SE.denom[1:(z-1)])
				non.censored <- length(which(!is.na(apply(xS,1,sum))))
			}
		}

		# Confidence limits
		cil.S.t <- S.t - (qt(0.975, df=M) * SE.t)
		ciu.S.t <- S.t + (qt(0.975, df=M) * SE.t)

		# output list
		out <- list(
			x = x, S=S, w=w,
			S.t = S.t, S.f = S.f, SE.it = SE.it, SE.t = SE.t, CI = rbind(CIL = cil.S.t, CIU = ciu.S.t),
			n.bar = n.bar, broods = counts, valid.broods = valid.counts,
			M = M, non.censored = non.censored, undefined = undefined)
		
		class(out) <- "flint"
		out
}

################################################################################
# x is a data frame with number of individuals in each brood
# g is a grouping variable (factor)
surv.prob.n <- function(x, g) {
	x.split <- split(x, g)
	lapply(x.split, surv.prob)
}

# Calculate test statistic D^2 for difference between two samples [survival functions]
d2.obs <- function(x, g) {
	if (length(levels(g)) > 2) stop("Only two groups can be compared!")
	surv.0 <- surv.prob.n(x=x, g=g)
	S.t.0 <- t(sapply(surv.0, "[[", "S.t"))
	SE.t.0 <- t(sapply(surv.0, "[[", "SE.t"))
	# Calculate test-statistic
	# z-score for each individual conditional survival time
	z0.t <- diff(apply(S.t.0, 2, rev)) / colSums(SE.t.0)
	# z.t <- diff(S.t.0) / colSums(SE.t.0)
	d2.obs <- sum(z0.t^2)
	list('S.t'=S.t.0, 'SE.t'=SE.t.0, 'z'=z0.t, 'd2'=d2.obs)
}
################################################################################


