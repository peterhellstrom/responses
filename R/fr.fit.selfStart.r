# Self-starting functions for estimation of non-linear functional responses
# Type IIIa & Type IIIb are not stable at all...

# Type II
#' @export
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

#' @export
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
#' @export
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

#' @export
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
#' @export
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

#' @export
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
#' @export
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

#' @export
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
#' @export
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

#' @export
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
#' @export
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
