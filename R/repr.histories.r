repr.histories <- function(x) {
	# Input x is a n-dimensional array created by table()
	x <- as.data.frame(x)
	x <- x[x$Freq > 0,]
	rownames(x) <- NULL
	x
}
