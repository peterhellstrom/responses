
########################################################################
################################################################################
# fResponse
# Functions for simulating and fitting functional responses
################################################################################
# Developed by Peter Hellstr?m, Department of Zoology, Stockholm University
# 2005-2013
# Last updated: November 22nd 2013
################################################################################
# NOTE: It is difficult to get nls to converge for the TypeIIIb function
# for some data sets. By using mle, I have realized that non-convergence for nls
# occurs when the mle-estimate for the exponent theta is less than one. That functional
# response is not sigmoid, but looks like a Type II response. So you better be careful
# when looking at AICc weights - if the Type IIIb response gets the highest weight
# one should also look at the values of the estimated parameter theta. It could be
# that the functional response is best described by a Type IIIb curve, but is not
# sigmoid in shape!
# It also seems that mle finds sensible estimates when nls fails. However,
# mle2 is slower than nls.



########################################################################
# Required packages for use of fResponse
# Check that required packages are installed

is.installed <- function(mypkg) mypkg %in% installed.packages()[,1]
install.missing <- function(mypkg, repos="http://ftp.acc.umu.se/mirror/CRAN/", dependencies=TRUE) {
	for (i in 1:length(mypkg)) if (is.installed(mypkg[i]) == FALSE) install.packages(pkgs=mypkg[i], lib=.Library, repos=repos, dependencies=dependencies)
}

pkgs <- c("bbmle","MASS","emdbook","plyr","nlstools","robustbase","nnet")
install.missing(pkgs)

# Load required packages
require(bbmle)
require(emdbook)
require(MASS)
require(plyr)
require(nlstools)
require(robustbase)
require(nnet)



########################################################################
# Stand alone version for calculating AIC and AICc values
# Use this function to compare AICc values
# AICc-functions are also available in the packages bbmle, MuMIn, AICcmodavg

# NOTE: beware that the attribute nobs differs between packages and functions
# nls & gnls return the actual number of observations
# gls return the number of observations - the number of model paramters (not including variance(s))!

aicc <- function(object, nobs=NULL) {
	# Get the number of datapoints and estimated parameters for different object classes
	# object.class <- class(object)[1]
	# cat(paste("class attribute of object:",object.class), "\n")
	if (!is.null(nobs)) {
		if(length(nobs) > 1) stop("nobs must be a scalar!")
		n <- nobs
	} else {
		n <- attr(logLik(object), "nobs")
	}
	k <- attr(logLik(object), "df")
	aic <- -2*logLik(object) + 2*k # AIC
	penalty <- (2*k*(k+1)) / (n-k-1) # penalty term
	if (n-k-1 == 0) cat("Warning: denominator of penalty term equals 0.\nAICc is infinite.\n")
	AICc <- aic + penalty
	# Create output vector
	out <- c('logLik'=logLik(object), 'k'=k, 'n'=n, 'AIC'=aic, 'AICc'=AICc, 'penalty'=penalty)
	out
}

# Extract AICc info and calculate statistics for many data sets
aicc.n <- function(object, sort=TRUE, round=FALSE, nobs=NULL, digits=5) {
	# Input is a list with objects (i.e. fitted models)
	n <- length(object)
	
	if (!is.null(nobs)) {
		tab <- data.frame(t(sapply(object,aicc,nobs)))
		tab$n <- nobs
	} else {
		tab <- data.frame(t(sapply(object,aicc)))
	}
	rownames(tab) <- names(object)
	index <- which(tab[,"AICc"] == min(tab[,"AICc"])) # Find model with smallest AICc
	
	deltai <- tab[,"AICc"] - min(tab[,"AICc"]) # AICc-delta values
	rel.like <- exp(-deltai / 2) # Relative likelihood
	wi <- rel.like / sum(rel.like) # Akaike weights
	ER <- wi[index] / wi # Information ratio
	ranking <- rank(deltai)
	
	# Create output (data frame)
	out <- data.frame('logLik' = tab$logLik, 'K' = tab$k, 'n'=tab$n, 'AIC' = tab$AIC, 'AICc' = tab$AICc,
		'delta' = deltai, 'weigths' = wi, 'ER' = ER, 'rank' = ranking)
	rownames(out) <- names(object)
	
	if (sort) out <- out[order(out[,"rank"]),]
	if (round) {
		round(out,digits)
	} else
	out
}



########################################################################
################################################################################
# Calculate number of killed prey (NP):
# n = number of individuals in each predator class (adults, juveniles etc.)
# Example for a raptor pair: c(adults, brood size), e.g. c(2,4)
# der = daily energy requirement in g
# pp = proportion of prey type in the diet (usually expressed in terms of biomass!)
# mmp = mean mass of prey type (in grams)
# T = Time, duration of study in days (defaults to 1)
# This function can use single numbers as inputs, or vectors (all vector must be of the same length)

NP <- function(n, der, pp, mmp, T=1) {
	inp <- list(der, pp, n)
	crit <- unique(sapply(inp, length))
	if (length(crit) > 1) stop("Input vectors der, pp & n must all be of equal length")
	T * (sum(n * der * pp) /  mmp)
}

################################################################################
# Daily energy requirement
# cast.rate = daily casting rate of pellets
# corr.factor = single value or vector with correction factor (corrects for digested prey, usually based on feeding experiments in captivity)
# mmp = mean mass in grams of prey tyep/species i 
# n = number of prey items of type/species i
# N = total number of pellets analyzed
# output = c("sum","ind"). If "sum" only the sum of all prey types is returned. If "ind", the contribution of each prey type i will be returned
der <- function(cast.rate=1.1, corr.factor=1, mmp, n, N, output=c("sum", "ind")) {
	output <- match.arg(output)
	if (output == "sum") {
		sum((cast.rate * corr.factor * mmp * n) / N)
	} else {
		(cast.rate * corr.factor * mmp * n) / N
	}
}

# Examples:
# der(mmp=c(31.50, 48.60, 31.30), n=c(10,15,14), N=20)
# Use different correction factors: input a vector of same length as n
# der(mmp=c(31.50, 48.60, 31.30), corr.factor=c(1.05,1.52,1.7), n=c(10,15,14), N=20)

################################################################################
# Predicted daily food requirements, from scaling theory
# After Nagy 1987 + Bozinovic & Medel 1988
# m = body weight (mass) of animal (in grams)
# eff = assimilation efficiency, a proportion (between 0 and 1)
# cal = caloric content of 1 g prey

# Default values: eff=0.769 (estimate for raptors), cal=6.65 (caloric content for 1 g of small mammal prey)
der.pred <- function(m=1000, eff=0.769, cal=6.65) {
	fmr <- 10.9 * m^0.64
	fmr / (eff * cal)
}

# Example:
# der.pred(m=800, eff=0.6)

################################################################################

# Heavier prey: fewer items are necessary to fulfil daily energy requirements.
# For heavier prey, asymptote of functional response will always be lower, so the relevant null hypotheses
# is not if asymptote differ between species, but rather if the asymptote differes more or less than expected!
# Expected difference is the weight ratio between prey types.

# Illustrated for grey-sided voles and lemmings, for the same proportion in diet:
# NP(pp=0.6089, der=140, mmp=31.50, n=4)
# NP(pp=0.6089, der=140, mmp=48.60, n=4)

# Input vectors
# NP(pp=c(0.5,0.5), der=c(110, 140), mmp=30, n=c(2,4))

# Example for female, male, nestling (females larger than males)
# NP(mmp=31.50, der=c(120,110,140), n=c(1,1,4), pp=c(0.6,0.6,0.75))

# Assume two predator classes/stages.
# If pooled proportion in diet in terms of biomass is p.hat=0.65, but proportions differ between classes.
# Assume prop1=0.8, then prop2 = 2*p.hat - prop1 = 2*0.65 - 0.8 = 0.5
# Compare estimates of NP, pooled vs. unpooled:
# NP(mmp=31.50, der=c(110,140), n=c(2,4), pp=c(0.65,0.65)) # Both classes have the same proportion = 16.10
 #NP(mmp=31.50, der=c(110,140), n=c(2,4), pp=c(0.5,0.8)) # Difference between prey classes = 17.71
# NP(mmp=31.50, der=c(110,140), n=c(2,4), pp=c(0.8,0.5)) # Switch proportion order = 14.47619
# Difference depends on der and n

################################################################################



########################################################################
extr.try <- function(x, fn) { 

	sets <- length(x)
	# rows is the number of succesful simulations
	# scc.res is a subset of res, containing only the succesful iterations:
	inds <- sapply(x, function(x) !inherits(x, "try-error"))
	scc.x <- x[inds]
	rows <- length(scc.x)
	n.success <- which(inds == TRUE)
	n.fail <- which(inds == FALSE)

	cat(paste("Length of input list: ", sets, "\n", sep=""))
	cat(paste("Number of successes: ", rows, "\n", sep=""))
	cat(paste("Failure: ", n.fail, "\n", sep=""))

	out <- sapply(scc.x, fn)
	out

}




########################################################################
# Convert between Holling's and Michaelis-Menten parameterization
convHollMM <- function (att,h) { 
	# att*x / (1 + att*h*x)
	a = 1/h
	b <- 1/(att*h)
	c(a=a, b=b)
}

convMMHoll <- function(a,b,theta=1) {
	# a*x^theta / (b^theta + x^theta)
	att = (a/b)^(1/theta)
	h = 1/a
	c(att=att, h=h)
}



########################################################################

# Update code to only include the multi.disc with minimization of the Box-Draper criterion!
# Procedure similar to current LLmulti.disc2
# LLmulti.disc only estimates response to one prey species
# LLmulti.disc2 estimates parameters to both species and is the preferred approach.

LLmulti.disc <- function(y, x1, x2, a1, a2, h1, h2) {
inds <- which(y<0) # do not allow negative values
if (length(inds)>0) {
  y <- y[-inds]
  x1 <- x1[-inds]
  x2 <- x2[-inds]
  } 

y.hat <- multi.disc(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2)$f1
  ssq <- sum((y.hat-y)^2)
  n <- length(y)
  sigma <- sqrt(ssq/n)
  return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

#####################
 
LLmulti.disc2 <- function(y1, y2, x1, x2, a1, a2, h1, h2) {
# I can't recall where this code came from! It is likely correct, but... 
  # Observed response
  y <- matrix(cbind(y1,y2),ncol=2,byrow=F)
  # Predicted response
  y.hat <- multi.disc(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2)
  mu <- matrix(cbind(y.hat$f1,y.hat$f2),ncol=2,byrow=F)
  # Residual matrix
  z <- t(y - mu)
  n <- length(mu)
  # Variance-covariance matrix of residuals
  Sigma <- (z %*% t(z))/n
  eS <- eigen(Sigma, sym = TRUE, EISPACK = TRUE)
  ev <- eS$values
  logdetS <- try(determinant(Sigma, logarithm = TRUE)$modulus)
  iS <- try(solve(Sigma))
  ssq <- diag(t(z) %*% iS %*% z)
  loglik = -(n * (log(2 * pi)) + logdetS + ssq)/2
  sum(loglik)

return(-sum(dmvnorm(x=y, mu=mu, Sigma=Sigma, log=TRUE)))
}

LLmulti.disc2m <- function(y1, y2, x1, x2, a1, a2, h1, h2, m1, m2) {
# I can't recall where this code came from! It is likely correct, but... 
  # Observed response
  y <- matrix(cbind(y1,y2),ncol=2,byrow=F)
  # Predicted response
  y.hat <- multi.disc.m(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2, m1=m1, m2=m2)
  mu <- matrix(cbind(y.hat$f1,y.hat$f2),ncol=2,byrow=F)
  # Residual matrix
  z <- t(y - mu)
  n <- length(mu)
  # Variance-covariance matrix of residuals
  Sigma <- (z %*% t(z))/n
  eS <- eigen(Sigma, sym = TRUE, EISPACK = TRUE)
  ev <- eS$values
  logdetS <- try(determinant(Sigma, logarithm = TRUE)$modulus)
  iS <- try(solve(Sigma))
  ssq <- diag(t(z) %*% iS %*% z)
  loglik = -(n * (log(2 * pi)) + logdetS + ssq)/2
  sum(loglik)

return(-sum(dmvnorm(x=y, mu=mu, Sigma=Sigma, log=TRUE)))
}

LLmulti.disc2bd <- function(y1, y2, x1, x2, a1, a2, h1, h2) {
  # Observed responses
  y <- matrix(c(y1,y2),ncol=2)
  # Predicted responses
  y.hat <- multi.disc(x1=x1, x2=x2, a1=a1, a2=a2, h1=h1, h2=h2)
  mu <- matrix(c(y.hat$f1,y.hat$f2),ncol=2)
  z <- y - mu
  # Box-Draper criterion (see chapter 4 in Bates & Watts 1988).
  # Minimizes the determinant of the crossproduct of the residual matrix.
  # Does not return a likelihood...but works still
  return(prod(svd(y - mu, nu=0, nv=0)$d)^2)
}

LLmulti.joly <- function(y, x1, x2, a.tot, h1, h2, alpha) {

y.hat <- multi.joly(x1=x1, x2=x2, a.tot=a.tot, h1=h1, h2=h2, alpha=alpha)
# Remove NA's
inds <- which(is.na(y.hat))

if (length(inds)>0) {
y <- y[-inds]
x1 <- x1[-inds]
x2 <- x2[-inds]
y.hat <- y.hat[-inds] }

ssq <- sum((y.hat-y)^2)
n <- length(y.hat)
sigma <- sqrt(ssq/n)
return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}



########################################################################
# Likelihood functions for functional responses
# All functions have normally distributed residuals
# Thus, the functions minimizes sums-of-squares and should give very similar results to nls()

# UPDATE 2013-11-14
# TODO: Check more about sigma...
# Previous versions did not include sigma as an estimated parameter
# I got this idea from the course at Grimsö with Tom Hobbs in 2008,
# but since then Tom has changed his course materials to include sigma as an estimated parameter.

# The original sigma statements are still in the code, but out-commented

# The new changes in these functions also have effects on fr.fit.r

# Null model, intercept only
LL.null <- function(x, a, sigma) {
	y.hat <- rep(mean(a),length(x))
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Type I, Linear response without intercept
LL.l <- function(y, x, a, sigma) {
	y.hat <- (a*x)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Michaelis-Menten, type II, IIIa & type IIIb
# Nested models, fix parameter theta=1 for Type II, theta=2 for Type IIIa
LL.mm <- function(y, x, a, b, theta, sigma) {
	#y.hat <- TypeIIIb(x, a, b, theta)
	y.hat <- fr.mm(x, a, b, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}

# Using weights
LL.mm.w <- function(y, x, weights, a, b, theta, sigma) {
	#y.hat <- TypeIIIb(x, a, b, theta)
	y.hat <- fr.mm(x, a, b, theta)
	# sigma <- sqrt(sum(weights*((y.hat-y)^2)) / length(x))
	return(-sum(dnorm(y,y.hat,sigma/sqrt(weights),log=TRUE)))
}

# Supply vector with parameters, not sure if I have ever used this approach
# Delete?
LL.mm.p <- function(p, x, y, sigma) {
	#if (is.null(names(p))) {
		#a <- p[1]
		#b <- p[2]
		#theta <- p[3]
	#} else {
		#if (any(c("a","b","theta") %in% names(p)==FALSE)) {
			#stop("Parameter names do not match")
		#} else {
			#a <- p[["a"]]
			#b <- p[["b"]]
			#theta <- p[["theta"]]
		#}
	#}
	a <- p[1]
	b <- p[2]
	theta <- p[3]
	
	y.hat <- fr.mm(x, a, b, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}


# Likelihood-function for two groups, type IIIb function
# Input variable: x,y,g [indexing variable for group membership]
LL.mm.2gr <- function(a1, a2, b1, b2, theta1, theta2, sigma) {
	g <- as.numeric(g)
	a <- c(a1,a2)[g]
	b <- c(b1,b2)[g]
	theta <- c(theta1, theta2)[g]
	y.hat <- TypeIIIb(x=x, a=a, b=b, theta=theta)
	# ssq <- sum((y.hat-y)^2)
	# n <- length(x)
	# sigma <- sqrt(ssq/n)
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE))) 
}

# Michaelis-Menten, type IVa & b
# Nested models, fix theta=1 for type IVa
LL.mm4 <- function(y, x, a, b, d, theta, sigma) {
	# y.hat <- TypeIVb(x, a, b, d, theta)
	y.hat <- fr.iv(x, a, b, d, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}	

# Holling's parameterization, type II, type IIIa & type IIIb
# Nested models, fix parameter theta to 1 for Type II, 2 for Type IIIa
LL.h <- function(y, x, a, h, theta, sigma) {
	# y.hat <- TypeIIIb.h(x, a, h, theta)
	y.hat <- fr.h(x, a, h, theta)
	# sigma <- sqrt(sum((y.hat-y)^2)/length(x))
	return(-sum(dnorm(y,y.hat,sigma,log=TRUE)))
}



########################################################################
# Multi-prey
# This code only estimates fr to one prey species,
# update to multi-response function (LLmulti.disc 2 or similar)

fr.fit.multi <- function(data) {

x1 <- data$x1
x2 <- data$x2
y1 <- data$y1
y2 <- data$y2

nlc <- nls.control(maxiter = 50000, tol=0.000001, minFactor=0.000001)

  inds <- which(x1<=0 | x2<=0 | y1<=0 | y2<=0)
  
  if (length(inds) > 0) {
  
  y1.new <- y1[-inds]
  y2.new <- y2[-inds]
  x1.new <- x1[-inds]
  x2.new <- x2[-inds]  

  mm.start.mle1 <- as.list(getInitial(y1.new ~ SSmicmen(x1.new, a, b),
              data=data.frame(x1.new,y1.new), control=nlc))
  start.mle1 <- list(a1=as.numeric((mm.start.mle1$a/mm.start.mle1$b)), 
              h1=as.numeric(1/mm.start.mle1$a))
                            
  mm.start.mle2 <- as.list(getInitial(y2.new ~ SSmicmen(x2.new, a, b),
              data=data.frame(x2.new,y2.new), control=nlc))
  start.mle2 <- list(a2=as.numeric((mm.start.mle2$a/mm.start.mle2$b)), 
              h2=as.numeric(1/mm.start.mle2$a))
              
  fm.multi <- mle2(minuslogl=LLmulti.disc2bd, start=
              list(a1=start.mle1$a1, a2=start.mle2$a2, h1=start.mle1$h1, h2=start.mle2$h2), 
              data=list(x1=x1,x2=x2,y1=y1,y2=y2), control=list(maxit=50000))
  }
  
  if (length(inds) == 0) {
  
  mm.start.mle1 <- as.list(getInitial(y1 ~ SSmicmen(x1, a, b),
              data=data.frame(x1,y1), control=nlc))
  start.mle1 <- list(a1=as.numeric((mm.start.mle1$a/mm.start.mle1$b)), 
              h1=as.numeric(1/mm.start.mle1$a))
                            
  mm.start.mle2 <- as.list(getInitial(y2 ~ SSmicmen(x2, a, b),
              data=data.frame(x2,y2), control=nlc))
  start.mle2 <- list(a2=as.numeric((mm.start.mle2$a/mm.start.mle2$b)), 
              h2=as.numeric(1/mm.start.mle2$a))
              
  fm.multi <- mle2(minuslogl=LLmulti.disc2bd, start=
              list(a1=start.mle1$a1, a2=start.mle2$a2, h1=start.mle1$h1, h2=start.mle2$h2), 
              data=list(x1=x1,x2=x2,y1=y1,y2=y2), control=list(maxit=5000))
  }
  
coefs <- rbind(
coef(fm.multi)
)

# colnames(coefs) <- c("a1","a2","h1","h2")
rownames(coefs) <- c("MultiDisc")

# Calculate AICc and AICc-weights:
# Function that calculates corrected AIC:
AICc <- function(object) {
k <- length(coef(object)) + 1
n <- length(x1)
AICc <- as.numeric(-2*logLik(object) + 2*k) + (2*k*(k+1)/(n-k-1))
out <- c(logLik(object),k,AICc)
names(out) <- c('LogL','K','AICc')
out
}

AICc.table <- c(
as.vector(AICc(fm.multi)[3])
)
names(AICc.table) <- c("MultiDisc")

deltai <- AICc.table - min(AICc.table)
rel.like <- exp(-deltai/2)
wi <- rel.like / sum(rel.like)
names(wi) <- c("MultiDisc")

fit.objects <- list(fm.multi=fm.multi)

out <- list(x1=x1, x2=x2, y=y1, start.mle1=start.mle1, start.mle2=start.mle2, 
coefs=coefs, AICc=AICc.table, AkaikeWeights=wi, fit.objects=fit.objects)
out

}



########################################################################
# Fit many functional response datasets
# Wrapper function for fr.fit

# Must be updated to include start values in fr.fit!!!
# This function DOES NOT WORK AT ALL yet due to changes in fr.fit

fr.fit.n <- function(data,method,eq) { 

	nsim <- length(data)

	# Repeat the curve-fitting "nsim" number of times set:
	# Store the output in a list called res:
	# Output: x,y, AICc and Akaike weights, stored in res - a huge list. 
	res <- lapply(1:nsim, function(i) try(fr.fit(data[[i]],method,eq),silent=TRUE))

	# Some generated data sets can not be fitted, so keep track of those:
	# rows is the number of succesful simulations
	# scc.res is a subset of res, containing only the succesful iterations:
	rows <- length(res[sapply(res, function(x) !inherits(x, "try-error"))])
	scc.res <- res[sapply(res, function(x) !inherits(x, "try-error"))]

	n.success <- which(sapply(res, function(x) !inherits(x, "try-error"))==TRUE)
	n.fail <- which(sapply(res, function(x) !inherits(x, "try-error"))==FALSE)

	# Extract Akaike weights
	Weights <- t(sapply(1:rows, function(i) scc.res[[i]]$AkaikeWeights))

	# FINAL output
	out <- list(
		out=scc.res,
		AICc.weights = Weights,
		iterations = nsim,
		succesful.iterations = rows,
		n.success=n.success,
		n.fail=n.fail
	)
	return(out)
}



########################################################################
# After the simulation: 
# Alternative power-analysis, check the number of times each model was selected
# as "best", using AICc as criterion:
# data is an object fitted with fr.fit or fr.fit.n

fr.fit.out <- function (data){
	# A function that finds the position of the maximum value in a vector
	max.AICc <- function(x) which.max(x)
	# Calculate in which column the maximum value is found, and make a table:
	tab1 <- table(apply(data$AICc.weights,1,max.AICc))
	# Make a better table
	inds <- as.numeric(names(tab1))
	tab2 <- numeric(5)
	names(tab2) <- c("Type0","TypeI","TypeII","TypeIIIa","TypeIIIb")
	tab2[inds] <- tab1
	# Calculate "Power", percentage of simulations that each model was selected
	pow <- round(tab2/sum(tab2),2)
	# Output
	list(n=sum(tab2),table=tab2,power=pow)
}



########################################################################
# Add a chosen type of functional response to an existing plot:
# (works only for objects fitted with nls, predict is not available for mle2.)
# Object = stored object, type = functional response type

plot.fResponse <- function(object, type=c("Type0","TypeI","TypeII","TypeIIIa","TypeIIIb","all"),
	main=NULL, xlab=NULL, ylab=NULL, lcols=c(1,1,2,3,4), ltys=c(3,1,2,3,4), lwd=2, dev.new=FALSE, ...) {
	
	if (length(lcols) < 5) lcols <- rep(lcols[1],5)
	if (length(ltys) < 5) ltys <- rep(ltys[1],5)
	coefs <- object$coefs
	
	if (missing(type)) {
		type <- "all"
	} else {
		type <- match.arg(type)
	}
	
	if (missing(main)) {
		strMain <- "Functional response"
	} else {
		strMain <- main
	}
	
	if (missing(xlab)) {
		strxlab <- "Prey density"
	} else {
		strxlab <- xlab
	}
	
	if (missing(ylab)) {
		strylab <- "Consumption rate"
	} else {
		strylab <- ylab
	}
	
	plot(object$x, object$y,
		bty="l", font.lab=2,
		xlab=strxlab, ylab=strylab, main=strMain,...)

	if (object$eq == "M-M") {
		if (type=="Type0") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			mtext(type)
		} else if (type=="TypeI") {
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			mtext(type)
		} else if (type=="TypeII") {
			with(as.list(coefs[3,]), curve(TypeII(x, a=a, b=b), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIa") {
			with(as.list(coefs[4,]), curve(TypeIIIa(x, a=a, b=b), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIb") {
			with(as.list(coefs[5,]), curve(TypeIIIb(x, a=a, b=b, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			mtext(type)
		} else if (type=="all") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			with(as.list(coefs[3,]), curve(TypeII(x, a=a, b=b), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			with(as.list(coefs[4,]), curve(TypeIIIa(x, a=a, b=b), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			with(as.list(coefs[5,]), curve(TypeIIIb(x, a=a, b=b, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			legend("bottomright",c("Type0","Type I","Type II","Type IIIa", "Type IIIb"), col=lcols, lty=ltys, lwd=lwd, bty="n")
			mtext("All fitted types")
		}
	} else if (object$eq == "Holling") {
		if (type=="Type0") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			mtext(type)
		} else if (type=="TypeI") {
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			mtext(type)
		} else if (type=="TypeII") {
			with(as.list(coefs[3,]), curve(TypeII.h(x, a=a, h=h), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIa") {
			with(as.list(coefs[4,]), curve(TypeIIIa.h(x, a=a, h=h), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIb") {
			with(as.list(coefs[5,]), curve(TypeIIIb.h(x, a=a, h=h, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			mtext(type)
		} else if (type=="all") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			with(as.list(coefs[3,]), curve(TypeII.h(x, a=a, h=h), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			with(as.list(coefs[4,]), curve(TypeIIIa.h(x, a=a, h=h), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			with(as.list(coefs[5,]), curve(TypeIIIb.h(x, a=a, h=h, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			legend("bottomright",c("Type 0", "Type I","Type II","Type IIIa", "Type IIIb"), col=lcols, lty=ltys, lwd=lwd, bty="n")
			mtext("All fitted types")
		}
	}
}



########################################################################
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



########################################################################
# Self-starting functions for estimation of non-linear functional responses
# Type IIIa & Type IIIb are not stable at all...

# Type II
SStypeII <- deriv(~ (a*x) / (b + x), c("a","b"), function(x,a,b) {})
attr(SStypeII,"pnames") <- c("a","b","sigma")
attr(SStypeII,"class") <- "selfStart"
attr(SStypeII,"initial") <- 
function (mCall, data, LHS) 
{
	# data = data frame with x and y data
	# LHS = left hand side = response variable = y
	# mCall = x
	# xy = sorted xy-data, sorted after x
    xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
	# Check number of data points
    if (nrow(xy) < 3) {
        stop("too few distinct input values to fit a Michaelis-Menten model")
    }
	# Lineweaver-Burk estimates
    pars <- as.vector(coef(lm(1/y ~ I(1/x), data = subset(xy, x!=0 & y!=0))))
	# Partially linear algorithm
	m1 <- nls(y ~ x/(b + x), data = xy, start = list(b = abs(pars[2L]/pars[1L])), algorithm = "plinear")
    # Extract coefficients
	# Calculate nls-estimate of sigma
	sigma <- as.numeric(sqrt(deviance(m1)/df.residual(m1)))
	pars <- as.vector(c(coef(m1),sigma))
	value <- c(pars[2L], pars[1L], pars[3L])
    names(value) <- c("a", "b", "sigma")
    value
}

SStypeII.nls <- deriv(~ (a*x) / (b + x), c("a","b"), function(x,a,b) {})
attr(SStypeII.nls,"pnames") <- c("a","b","sigma")
attr(SStypeII.nls,"class") <- "selfStart"
attr(SStypeII.nls,"initial") <- 
function (mCall, data, LHS) 
{
	# data = data frame with x and y data
	# LHS = left hand side = response variable = y
	# mCall = x
	# xy = sorted xy-data, sorted after x
    xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
	# Check number of data points
    if (nrow(xy) < 3) {
        stop("too few distinct input values to fit a Michaelis-Menten model")
    }
	# Lineweaver-Burk estimates
    pars <- as.vector(coef(lm(1/y ~ I(1/x), data = subset(xy, x!=0 & y!=0))))
	# Partially linear algorithm
	m1 <- nls(y ~ x/(b + x), data = xy, start = list(b = abs(pars[2L]/pars[1L])), algorithm = "plinear")
    # Extract coefficients
	pars <- as.vector(coef(m1))
	value <- c(pars[2L], pars[1L])
    names(value) <- c("a", "b")
    value
}

# Type IIIa
SStypeIIIa <- deriv(~ (a*x^2) / (b^2 + x^2), c("a","b"), function(x,a,b) {})
attr(SStypeIIIa,"pnames") <- c("a","b","sigma")
attr(SStypeIIIa,"class") <- "selfStart"
attr(SStypeIIIa,"initial") <- 
function (mCall, data, LHS) 
{
	# data = data frame with x and y data
	# LHS = left hand side = response variable = y
	# mCall = x
	# xy = sorted xy-data, sorted after x
    xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
	# Check number of data points
    if (nrow(xy) < 3) {
        stop("too few distinct input values to fit a model")
    }
	pars <- as.vector(coef(lm(1/y ~ I(1/x), data = subset(xy, x!=0 & y!=0))))
	p <- list(b = abs(pars[2L]/pars[1L]))
	m1 <- nls(y ~ x^2/(b^2 + x^2), data = xy, start = p, algorithm="plinear")
	# Extract coefficients
	# Calculate nls-estimate of sigma
	sigma <- as.numeric(sqrt(deviance(m1)/df.residual(m1)))
	pars <- as.vector(c(coef(m1),sigma))
	value <- c(pars[2L], pars[1L], pars[3L])
    names(value) <- c("a", "b", "sigma")
    value
}

SStypeIIIa.nls <- deriv(~ (a*x^2) / (b^2 + x^2), c("a","b"), function(x,a,b) {})
attr(SStypeIIIa.nls,"pnames") <- c("a","b")
attr(SStypeIIIa.nls,"class") <- "selfStart"
attr(SStypeIIIa.nls,"initial") <- 
function (mCall, data, LHS) 
{
	# data = data frame with x and y data
	# LHS = left hand side = response variable = y
	# mCall = x
	# xy = sorted xy-data, sorted after x
    xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
	# Check number of data points
    if (nrow(xy) < 3) {
        stop("too few distinct input values to fit a model")
    }
	pars <- as.vector(coef(lm(1/y ~ I(1/x), data = subset(xy, x!=0 & y!=0))))
	p <- list(b = abs(pars[2L]/pars[1L]))
	m1 <- nls(y ~ x^2/(b^2 + x^2), data = xy, start = p, algorithm="plinear")
	# Extract coefficients
	pars <- as.vector(coef(m1))
	value <- c(pars[2L], pars[1L])
    names(value) <- c("a", "b")
    value
}

# Type IIIb
SStypeIIIb <- deriv(~ (a*x^theta) / (b^theta + x^theta), c("a","b","theta"), function(x,a,b,theta) {})
attr(SStypeIIIb,"pnames") <- c("a","b","theta","sigma")
attr(SStypeIIIb,"class") <- "selfStart"
attr(SStypeIIIb,"initial") <- 
function (mCall, data, LHS) 
{
	# data = data frame with x and y data
	# LHS = left hand side = response variable = y
	# mCall = x
	# xy = sorted xy-data, sorted after x
	# xy <- data.frame(sortedXyData(x, y, data))
    xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
	# Check number of data points
    if (nrow(xy) < 3) {
        stop("too few distinct input values to fit a model")
    }
	nlc <- nls.control(maxiter = 5000, tol=0.000001, minFactor=0.000001)
	# try to use type II-Lineweaver-Burk estimate for b, and default value of theta=2
	# Can't figure out a robust way to estimate theta without logs etc.
	pars <- as.vector(coef(lm(1/y ~ I(1/x), data = subset(xy, x!=0 & y!=0))))
	p <- list(b = abs(pars[2L]/pars[1L]), theta = 2)
	m1 <- nls(y ~ (x^theta)/(b^theta + x^theta), data = xy, start = p, algorithm="plinear", control=nlc)
	# Extract coefficients
	# Calculate nls-estimate of sigma
	sigma <- as.numeric(sqrt(deviance(m1)/df.residual(m1)))
	pars <- as.vector(c(coef(m1),sigma))
	value <- c(pars[3L], pars[1L], pars[2L], pars[4L])
    names(value) <- c("a", "b", "theta", "sigma")
    value
}

SStypeIIIb.nls <- deriv(~ (a*x^theta) / (b^theta + x^theta), c("a","b","theta"), function(x,a,b,theta) {})
attr(SStypeIIIb.nls,"pnames") <- c("a","b","theta")
attr(SStypeIIIb.nls,"class") <- "selfStart"
attr(SStypeIIIb.nls,"initial") <- 
function (mCall, data, LHS) 
{
	# data = data frame with x and y data
	# LHS = left hand side = response variable = y
	# mCall = x
	# xy = sorted xy-data, sorted after x
	# xy <- data.frame(sortedXyData(x, y, data))
    xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
	# Check number of data points
    if (nrow(xy) < 3) {
        stop("too few distinct input values to fit a model")
    }
	nlc <- nls.control(maxiter = 5000, tol=0.000001, minFactor=0.000001)
	# try to use type II-Lineweaver-Burk estimate for b, and default value of theta=2
	# Can't figure out a robust way to estimate theta without logs etc.
	pars <- as.vector(coef(lm(1/y ~ I(1/x), data = subset(xy, x!=0 & y!=0))))
	p <- list(b = abs(pars[2L]/pars[1L]), theta = 2)
	m1 <- nls(y ~ (x^theta)/(b^theta + x^theta), data = xy, start = p, algorithm="plinear", control=nlc)
	# Extract coefficients
	pars <- as.vector(coef(m1))
	value <- c(pars[3L], pars[1L], pars[2L])
    names(value) <- c("a", "b", "theta")
    value
}

# code chunk previously used to estimate start values for sigmoid responses

	#a.hat <- as.numeric(ceiling(max(xy$y)))
	#z <- xy[["y"]]
	#if (min(z) <= 0) z <- abs(z)
	#xy[["z"]] <- log(z/(a.hat - z))
	
	#m0 <- lm(z ~ log(x), data=xy)
	#theta.hat <- as.numeric(coef(m0)[2]) # slope
	#b.hat <- as.numeric(exp(-coef(m0)[1] / theta.hat))
	#p <- list(b = b.hat, theta=theta.hat)

# Wrapper function to generate start values easier

# Must also include estimate of sigma for mle2-functions
# getInitial(y ~ SStypeIIIb(x,a,b,theta), data=xy.list) part is wrong, this does not work...
getStart.mm <- function(x, y, method=c("mle","nls")) {
	
	method = match.arg(method)
	xy.list <- list(x=x,y=y)

	if (method == "mle") {
		ini.0 <- c('a'=mean(y), 'sigma'=sd(y))
		m.I.start <- lm(y~x-1)
		ini.I <- c('a'=as.vector(coef(m.I.start)), 'sigma'=as.numeric(sqrt(deviance(m.I.start)/df.residual(m.I.start))))
	
		out <- list(
			'eq' = 'M-M',
			'x' = x,
			'y' = y,
			'ini.0' = ini.0,
			'ini.I' = ini.I,
			'ini.II' = getInitial(y ~ SStypeII(x,a,b), data=xy.list),
			'ini.IIIa' = getInitial(y ~ SStypeIIIa(x,a,b), data=xy.list),
			'ini.IIIb' = getInitial(y ~ SStypeIIIb(x,a,b,theta), data=xy.list)
		)
	} else if (method == "nls") {
		ini.0 <- c('a'=mean(y))
		m.I.start <- lm(y~x-1)
		ini.I <- c('a'=as.vector(coef(m.I.start)))
	
		out <- list(
			'eq' = 'M-M',
			'x' = x,
			'y' = y,
			'ini.0' = ini.0,
			'ini.I' = ini.I,
			'ini.II' = getInitial(y ~ SStypeII.nls(x,a,b), data=xy.list),
			'ini.IIIa' = getInitial(y ~ SStypeIIIa.nls(x,a,b), data=xy.list),
			'ini.IIIb' = getInitial(y ~ SStypeIIIb.nls(x,a,b,theta), data=xy.list)
		)
	
	}
	
	class(out) <- "getStart"
	out
}

getStart.mm2 <- function(x, y, method=c("mle","nls"), theta=2) {
	# Use the estimates from type IIIa for type IIIb as well!
	
	method <- match.arg(method)
	xy.list <- list(x=x,y=y)

	if (method == "mle") {
	
		ini.0 <- c('a'=mean(y), 'sigma'=sd(y))
		m.I.start <- lm(y~x-1)
		ini.I <- c('a'=as.vector(coef(m.I.start)), 'sigma'=as.numeric(sqrt(deviance(m.I.start)/df.residual(m.I.start))))
		ini.IIIb <- c(getInitial(y ~ SStypeIIIa(x,a,b), data=xy.list), 'theta'=theta)
		ini.IIIb <- ini.IIIb[c(1,2,4,3)]
	
		out <- list(
			'eq' = 'M-M',
			'x' = x,
			'y' = y,
			'ini.0' = ini.0,
			'ini.I' = ini.I,
			'ini.II' = getInitial(y ~ SStypeII(x,a,b), data=xy.list),
			'ini.IIIa' = getInitial(y ~ SStypeIIIa(x,a,b), data=xy.list),
			'ini.IIIb' = ini.IIIb
		)
	} else if (method == "nls") {
		ini.0 <- c('a'=mean(y))
		m.I.start <- lm(y~x-1)
		ini.I <- c('a'=as.vector(coef(m.I.start)))
		ini.IIIb <- c(getInitial(y ~ SStypeIIIa.nls(x,a,b), data=xy.list), 'theta'=theta)
		ini.IIIb <- ini.IIIb[c(1,2,4,3)]
	
		out <- list(
			'eq' = 'M-M',
			'x' = x,
			'y' = y,
			'ini.0' = ini.0,
			'ini.I' = ini.I,
			'ini.II' = getInitial(y ~ SStypeII.nls(x,a,b), data=xy.list),
			'ini.IIIa' = getInitial(y ~ SStypeIIIa.nls(x,a,b), data=xy.list),
			'ini.IIIb' = ini.IIIb
		)
	
	}
	
	class(out) <- "getStart"
	out
}

# Add function that calculate Holling-type start parameters from Michaelis-Menten parameterization:
# Check how these values really should be calculated
# Something is wrong - I can't use these estimates even if I simulate from a well behaved disk-equation...
getStart.h <- function(x, y, method=c("mle","nls"), theta=2) {
	method <- match.arg(method)
	xy.list = list(x,y)
	
	if (method == "mle") {
		mm <- getStart.mm2(x, y, theta, method="mle")
		ini.II = with(as.list(mm$ini.II), c('a'=a/b, 'h'=1/a, 'sigma'=sigma))
		ini.IIIa <- with(as.list(mm$ini.IIIa), c('a'=a/b^2, 'h'=1/a, 'sigma'=sigma))
		ini.IIIb <- with(as.list(mm$ini.IIIb), c('a'=2 * (a/b^theta), 'h'=1/a, 'theta'=theta, 'sigma'=sigma))
		out <- list('eq'='Holling', 'x' = x,'y' = y, 'ini.0' = mm$ini.0, 'ini.I' = mm$ini.I, 'ini.II'=ini.II, 'ini.IIIa'=ini.IIIa, 'ini.IIIb'=ini.IIIb)
	} else if (method == "nls") {
		mm <- getStart.mm2(x, y, theta, method="nls")
		ini.II = with(as.list(mm$ini.II), c('a'=a/b, 'h'=1/a))
		ini.IIIa <- with(as.list(mm$ini.IIIa), c('a'=a/b^2, 'h'=1/a))
		ini.IIIb <- with(as.list(mm$ini.IIIb), c('a'=2 * (a/b^theta), 'h'=1/a, 'theta'=theta))
		out <- list('eq'='Holling', 'x' = x,'y' = y, 'ini.0' = mm$ini.0, 'ini.I' = mm$ini.I, 'ini.II'=ini.II, 'ini.IIIa'=ini.IIIa, 'ini.IIIb'=ini.IIIb)
	
	}
	class(out) <- "getStart"
	out
}

getStart.method <- function(x, y, plot=TRUE) {

	m1 <- lm(x/y ~ x)
	m2 <- lm(1/y ~ I(1/x))
	
	bolker <- c('att'=as.numeric((coef(m1)[1])^-1), 'h'=as.numeric(coef(m1)[2]))
	bolker.mm <- convHollMM(att=as.numeric(bolker["att"]), h=as.numeric(bolker["h"]))
	
	out <- list(
		'bolker' = bolker,
		'bolker.mm' = bolker.mm,
		'lineweaver-burk' = c('a' = as.numeric(1/abs(coef(m2)[1])), 'b'= as.numeric(abs(coef(m2)[2]/coef(m2)[1])))
	)
	
	if (plot) {
		par(mfrow=c(2,2))
			plot(x,x/y)
			plot(1/x, 1/y)
			plot(x,x/y)
				abline(m1,col=2)
			plot(1/x,1/y)
				abline(m2,col=2)
			par(mfrow=c(1,1))
	}
	
	out
}


# Function that plots initial estimates for a given (single) dataset
plot.getStart <- function(object, ...) {
	
	eq <- object$eq
	x <- object$x
	y <- object$y
	plot(x,y,...)
	
	if (eq=="M-M") {
		with(as.list(object$ini.II), curve(TypeII(x, a=a, b=b), add=TRUE, col=2,lty=2))
		with(as.list(object$ini.IIIa), curve(TypeIIIa(x, a=a, b=b), add=TRUE, col=3, lty=3))
		with(as.list(object$ini.IIIb), curve(TypeIIIb(x, a=a, b=b,theta=theta), add=TRUE, col=4, lty=4))
	} else if (eq == "Holling") {
		with(as.list(object$ini.II), curve(TypeII.h(x, a=a, h=h), add=TRUE, col=2,lty=2))
		with(as.list(object$ini.IIIa), curve(TypeIIIa.h(x, a=a, h=h), add=TRUE, col=3, lty=3))
		with(as.list(object$ini.IIIb), curve(TypeIIIb.h(x, a=a, h=h, theta=theta), add=TRUE, col=4, lty=4))
	}
	legend("bottomright",c("Type II","Type IIIa", "Type IIIb"), col=2:4, lty=2:4, bty="n")
}



########################################################################
# Set up functions for the different response types:
# (Generating functions)

# SINGLE PREY
# Michaelis-Menten type
Type0 <- function(x,a) rep(mean(a),length(x))
TypeI <- function(x,a) (a*x)
TypeI.intercept <- function(x,a,b) (a*x + b)

TypeII <- deriv(~ (a*x) / (b + x), c("a","b"), function(x,a,b) {})
TypeII.shifted <- deriv(~ ((a*(x-c)) / (b + (x-c))) + d, c("a","b","c","d"), function(x,a,b,c,d) {})
TypeIIIa <- deriv(~ a*x^2 / (b^2 + x^2), c("a","b"), function(x,a,b) {})
TypeIIIb <- deriv(~ a*x^theta / (b^theta + x^theta), c("a","b","theta"), function(x,a,b,theta) {})

TypeIVa <- deriv(~ a*x^2 / (b+d*x + x^2), c("a","b","d"), function(x,a,b,d) {})
TypeIVb <- deriv(~ a*x^theta / (b+d*x + x^theta), c("a","b","d","theta"), function(x,a,b,d,theta) {})

# Original Holling parameterization (a & h)
TypeII.h <- deriv(~ (a*x / (1 + a*h*x)), c("a","h"), function(x,a,h) {})
TypeIIIa.h <- deriv(~ (a*x^2) / (1 + a*h*x^2), c("a","h"), function(x,a,h) {})
TypeIIIb.h <- deriv(~ (a*x^theta) / (1 + a*h*x^theta), c("a","h","theta"), function(x,a,h,theta) {})

# Functions to use together with mle2 (the above functions do not always work with mle2, probably because the deriv/gradient part).
fr.mm <- function(x,a,b,theta) a*x^theta / (b^theta + x^theta)
fr.iv <- function(x,a,b,d,theta) a*x^theta / (b+d*x + x^theta)
fr.h <- function(x,a,h,theta) (a*x^theta) / (1 + a*h*x^theta)

# NOT YET IMPLEMENTED IN THE UPDATE fResponse.r
# MULTI-PREY

# Should be updated to:
# 1) Allow for more than two species
# 2) Use input vectors, e.g. x=c(x1,x2), a=c(a1,a2), h=c(h1,h2), theta=c(theta1,theta2)
# 3) The denominator is common in all cases, calculate only once with matrix multiplication.
# This is done for TypeII.h.n, not yet for the other two...

multi.disc <- function(x1,x2,a1,a2,h1,h2) {
	f1 <- a1*x1 / (1 + a1*h1*x1 + a2*h2*x2)
	f2 <- a2*x2 / (1 + a1*h1*x1 + a2*h2*x2)
	list(f1=f1,f2=f2)
}

TypeII.h.n <- function(x, a, h, output=c("list","matrix")) {
	output <- match.arg(output)
	# x = matrix with n column
	# a = vector with attack rates of length n
	# h = vector with handling times of length n
	# output = c("list", "matrix"): output in list or matrix format
	
	if (ncol(x) != length(a) | ncol(x) != length(h)) stop("Wrong dimensions")
	
	num <- sweep(x,2,a,'*')
	denom <- 1 + (x %*% (a*h))
	f <- sweep(num,1,denom,'/')
	labs <- paste("f",1:ncol(x),sep="")
	
	if (output == "matrix") {
		colnames(f) <- labs
	} else if (output == "list") {
		f <- split(f, rep(1:ncol(f), each = nrow(f)))
		names(f) <- labs
	}
	
	f
}

multi.disc.m <- function(x1,x2,a1,a2,h1,h2,m1,m2) {
	f1 <- a1*x1^m1 / (1 + a1*h1*x1^m1 + a2*h2*x2^m2)
	f2 <- a2*x2^m2 / (1 + a1*h1*x1^m1 + a2*h2*x2^m2)
	list(f1=f1,f2=f2)
}

multi.joly <- function(x1,x2,a.tot,h1,h2,alpha) {
	f1 <- a.tot*alpha*x1 / (1 + a.tot*alpha*h1*x1 + a.tot*(1-alpha)*h2*x2)
	f1
}




########################################################################
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



########################################################################
fr.sim.data.multi <- function(parms,lims,sigma,n,sets,cor.method,...) {

# Simulate two-prey systems, mainly for testing estimation functions
# Works like:
# 1) create prey densities, draw from uniform distribution, or from a multivariate normal
# (option cor.method). Drawing from a multivariate normal distribution allows for correlation
# between prey numbers, but a variance-covariance matrix has to be defined for use.
# 2) Calculate the deterministic (y.hat) value, based on the values in the parms argument
# 3) Calculate the stochastic value (observation error) based on y.hat ("mean"), by drawing from
# a normal distribution with sd=sigma argument.
# cova = covariance matrix???
 
out <- vector("list",sets)

for (i in 1:sets) {
  
  if (cor.method=="null") {
  
  x1 <- runif(n=n, min=lims$min.x1, max=lims$max.x1)
  x2 <- runif(n=n, min=lims$min.x2, max=lims$max.x2)
   
  }
  
  if (cor.method=="correlated") {
  
  cova <- cova

  x <- mvrnorm(n, mu = c(cor.parms$x1,cor.parms$x2), Sigma = cova)

  inds <- which(x[1:n,]<0)
  x[inds] <- NA
  
  x1 <- x[,1]
  x2 <- x[,2]
  }
  
  y.hats <- multi.disc(x1=x1,x2=x2,a1=parms$a1,a2=parms$a2,h1=parms$h1,h2=parms$h2)
   
  y1.hat <- y.hats$f1
  y1 <- abs(rnorm(x1,mean=y1.hat,sd=sigma))
  # inds1 <- which(y1<0)
  # y1[inds1] <- NA
    
  y2.hat <- y.hats$f2
  y2 <- abs(rnorm(n,mean=y2.hat,sd=sigma))
  # inds2 <- which(y2<0)
  # y2[inds2] <- NA
  
  out[[i]] <- list(
  input = data.frame(a1=parms$a1, a2=parms$a2, h1=parms$h1, h2=parms$h2, min.x1=lims$min.x1, max.x1=lims$max.x1, 
          min.x2=lims$min.x2, max.x2=lims$max.x2, sigma=sigma, n=n, sets=sets),
  x1=x1,
  y1.hat=y1.hat,
  y1=y1,
  x2=x2,
  y2.hat=y2.hat,
  y2=y2)
  
  }
out
}



########################################################################
# Function for simulating data sets:
# Generate prey density (x), calculate deterministic y, and error normally distributed
# error to y.det. Parameterization is Michaelis-Menten.

# Inputs: 
# parms = functional response parameters a,b,c (must be a list)
# Functional response parameters:
# a = asymptote (miximum kill rate)
# b = half-saturation constant
# c = exponent controlling shape of curve
# lims = min and max x-values (must be called min.x and max.x, list)
# sigma = error (rse)
# n = number of data points in each set
# sets = generate "sets" number of data sets
# xmethod = how to generate x-values:
# full.x.range = draw n uniform values from min.x to max.x
# even.x.range = split x-range into 4 bins (determined by the quantiles of min.x and max.x)# , assign equal (or nearly equal if n isn't a multiple of 4) number of values in each bin.

fr.sim.data <- function(parms,lims,sigma,n,sets,xmethod) {

out <- vector("list",sets)

for (i in 1:sets) {

if (xmethod=="full.x.range") {   
  x <- runif(n=n, min=lims$min.x, max=lims$max.x) 
  }
  
if (xmethod=="even.x.range") {

  groups <- quantile(c(lims$min.x, lims$max.x))
  x.vec <- rep(1:4,times=ceiling(n/4))
  n.samp <- table(sample(x.vec,n,replace=FALSE))

  x <- c(
        runif(n.samp[1],min=groups[1],max=groups[2]),
        runif(n.samp[2],min=groups[2],max=groups[3]),
        runif(n.samp[3],min=groups[3],max=groups[4]),
        runif(n.samp[4],min=groups[4],max=groups[5])
        )
  }
  
  y.hat <- parms$a*x^parms$c / (parms$b^parms$c + x^parms$c)
  # Draw random normal variates and take absolute value
  # Not the best solution, should find a better solution
  # based on Poission or negbin
  y <- abs(rnorm(x,mean=y.hat,sd=sigma))
  
  out[[i]] <- list(
  input=data.frame(a=parms$a, b=parms$b, c=parms$c, min.x=lims$min.x, max.x=lims$max.x,       
        sigma=sigma, n=n, sets=sets),
  x=x,
  y.hat=y.hat,
  y=y)
  
  }
out
}



########################################################################
# A simple grid search function
# This is one way to find start-values for estimation of functional response parameters.
# NOTE: the nls2 package has more advanced grid search options for start values!

fr.start.grid <- function(x, y, type=c("TypeII","TypeIIIa","TypeIIIb"), length.out=10, extra=50) {

	type <- match.arg(type)

	# Generate grid
	# For
	# a (asymptote): search in interval mean(y) to max(y) + extra
	# b (half-saturation constant): 0 to mean(x)
	# theta: 0.5 to 5
	
	if (type == "TypeII") {
		grid <- expand.grid(a = seq(mean(y), max(y) + extra, length.out=length.out), b = seq(0, mean(x), length.out=length.out))
		ss <- function(p) sum((y - TypeII(x, p[1], p[2]))^2)
		idx  <- which.min(apply(grid, 1, ss))
		startval <- grid[idx,]
		fm <- nls(y ~ TypeII(x, a, b), start = startval)
	} else if (type == "TypeIIIa") {
		grid <- expand.grid(a = seq(mean(y), max(y) + extra, length.out=length.out), b = seq(0, mean(x), length.out=length.out))
		ss <- function(p) sum((y - TypeIIIa(x, p[1], p[2]))^2)
		idx  <- which.min(apply(grid, 1, ss))
		startval <- grid[idx,]
		fm <- nls(y ~ TypeIIIa(x, a, b), start = startval)
	} else if (type == "TypeIIIb") {
		grid <- expand.grid(a = seq(mean(y), max(y) + extra, length.out=length.out), b = seq(0, mean(x), length.out=length.out), theta = seq(0.5, 5, length.out=length.out))
		ss <- function(p) sum((y - TypeIIIb(x, p[1], p[2], p[3]))^2)
		idx  <- which.min(apply(grid, 1, ss))
		startval <- grid[idx,]
		fm <- nls(y ~ TypeIIIb(x, a, b, theta), start = startval)
	} else stop("Model not included!")

	dat <- list('start'=startval, 'nls'=fm, 'nlsSum'=summary(fm)) 
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
# Function that obtains likelihood-profiles and confidence intervals for nls- or mle2-objects
pci <- function(object, plot=FALSE, output=c("simple","full")) {
	
	if (class(object) == "lm") stop("object class lm is not implemented")
	else if (class(object) == "glm") stop("object class glm is not supported")

	#p <- profile(object)
	p <- object
	if (plot) plot(profile(object))
	output <- match.arg(output)
	
	if (output == "simple") {
		if (class(object)=="mle2") {
			out <- data.frame(
				'coef' = object@coef,
				'se' = summary(object)@coef[,2],
				'ci' = confint(p)
			)
		} else if (class(object)=="nls") {
			out <- data.frame(
				'coef' = coef(object),
				'se' = summary(object)$coef[,2],
				'ci' = confint(p)
				)
		}
	} else if (output == "full") {
		if (class(object)=="mle2") {
			out <- list(
				'coef' = object@coef,
				'vcov' = vcov(object),
				'cov2cor' = cov2cor(vcov(object)),
				'se' = summary(object)@coef[,2],
				'ci' = t(confint(p))
				)
		} else if (class(object)=="nls") {
			out <- list(
				'coef' = coef(object),
				'vcov' = vcov(object),
				'cov2cor' = cov2cor(vcov(object)),
				'se' = summary(object)$coef[,2],
				'ci' = t(confint(p))
			)
		}
	}
	
	out
}




########################################################################
# ANALYZE DENSITY-DEPENDENCE IN PREDATION RATE
# Fit polynomial regression of order 0-3 to the data.

# General function
pred.rate <- function(x, y, plot=TRUE, dev.new=TRUE) {
	
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
		
		
		if (dev.new) dev.new()
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



########################################################################
# Make a function that illustrates switching curves
# c = preference
# b = how preference changes with N1/N2, (b = 2, linear change)

switching <- function(c, b, type=c("proportions","ratios"), round=3, x.max=NULL, ...) {
	
	b.text = round(b,round)
	c.text = round(c,round)
	
	type <- match.arg(type)
	if (type == "proportions") {
		# Basic plot
		plot(x=c(0,1),y=c(0,1),type="n",
			xlab="Proportion available", ylab="Proportion in diet",
			font.lab=2,bty="l",las=1, ...)
		# Null hypothesis
		curve(switch.curve(x,b=1,c), n=1001, lty=5, lwd=2, col=2, add=TRUE)
		# Switching curve
		curve(switch.curve(x,b,c), n=1001, lty=1, lwd=2, col=1, add=TRUE)
		mtext(paste("Switching, c =",c.text,", b =",b.text))
		legend("bottomright",legend=c("Switching","Null case"),lty=c(1,5),lwd=c(2,2),col=c(1,2),bty="n",cex=1)
	}
	else if (type == "ratios") {
		# Basic plot
		y <- switch.curve.ratios(seq(0,x.max,0.1),b,c)
		plot(x=c(0,x.max),y=c(0,max(y)),type="n",
			xlab="Ratio available",ylab="Ratio in diet",
			font.lab=2,bty="l",las=1,ylim=c(0,max(y)), ...)
		# Null hypothesis
		curve(switch.curve.ratios(x,b=1,c), n=1001, lty=5, lwd=2, col=2, add=TRUE)
		# Switching curve
		curve(switch.curve.ratios(x,b,c), n=1001, lty=1, lwd=2, col=1, add=TRUE)
		mtext(paste("Switching, c =",c.text,", b =",b.text))
		legend("bottomright",legend=c("Switching","Null case"),lty=c(1,5),lwd=c(2,2),col=c(1,2),bty="n",cex=1)
	}
}

switch.curve <- function(x, b, c) c*x^b / ((1 - x)^b + c*x^b)
# x = proportion available
# c = preference
# b = how preference changes with N1/N2, (b = 2, linear change)
switch.curve.ratios <- function(x, b, c) c*x^b


