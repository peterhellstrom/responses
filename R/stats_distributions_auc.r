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

# Create vectors with x- and y-coordinates for area(s) under a given distribution.
# distr = name of distribution, in R's standard annotation for distributions, e.g. norm, pois, lnorm.
# para = a named list (names must correspond to arguments of distr) with parameter values.
# from, to = range in which to return polygon coordinates.
# by = increment

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
