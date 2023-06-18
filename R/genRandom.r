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
# .x <- genRandom(distr="norm", para=list(n=10000, mean=c(1,10,20), sd=c(2,2,2)), catgs=3, print=TRUE)
# hist(.x$x, col="lightgrey", breaks=60)

# A lot of the code is borrowed from the package fitdistrplus

genRandom <- function(distr, para, catgs=3, g=NULL, prob=NULL, print=FALSE, ...) {
	
	# Evaluate input arguments
	if(class(para) != "list") stop("'para' must be a named list")
	if (is.null(names(para))) stop("'para' must be a named list")
	if ((missing(distr) & !missing(para)) || (missing(distr) & !missing(para))) stop("'distr' and 'para' must defined")
	if (missing(catgs)) stop("'catgs', number of categories must be entered")
	if (!is.null(g)) if((class(g) %in% c("factor","integer")==FALSE)) stop("Grouping variable 'g' must be numeric")
	
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
	
	# Create grouping variable, if a grouping variable is not supplied by the user:
	if (is.null(g)) {
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
	if (print==TRUE) x <- data.frame(g,x)
	# Print
	x
}
