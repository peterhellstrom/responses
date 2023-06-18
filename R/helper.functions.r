# Function for writing list objects to text file(s)
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
