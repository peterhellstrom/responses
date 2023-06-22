# Different functional and numerical response types

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

# "Standard" type II & type III models
type2 <- deriv(~ (a*x) / (b + x), c("a","b"), function(x,a,b) {})
type3 <- deriv(~ a*x^2 / (b^2 + x^2), c("a","b"), function(x,a,b) {})
type3b <- deriv(~ a*x^theta / (b^theta + x^theta), c("a","b","theta"), function(x,a,b,theta) {})

# Threshold models
type2ta <- deriv(~ (a*(x-d)) / ((b-d) + (x-d)), c("a","b","d"), function(x,a,b,d) {})
type2tb <- deriv(~ (a*(x-d)) / (b + (x-d)), c("a","b","d"), function(x,a,b,d) {})

type3ta <- deriv(~ (a*(x-d)^theta) / ((b-d)^theta + (x-d)^theta), c("a","b","d","theta"), function (x,a,b,d,theta) {} )
type3tb <- deriv(~ (a*(x-d)^theta) / ((b)^theta + (x-d)^theta), c("a","b","d","theta"), function (x,a,b,d,theta) {} )

# Different parameterizations of the logistic function used by Henden et al.
logist <- deriv(~ a/(1 + exp(-(b + d*x))), c("a","b","d"), function(x,a,b,d) {})
logist2 <- deriv(~ a/(1 + exp((b + d*-x))), c("a","b","d"), function(x,a,b,d) {})
logist3 <- deriv(~ a/(1 + exp((b - d*x))), c("a","b","d"), function(x,a,b,d) {})

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

