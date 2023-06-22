#' @export
FreqTogen <- function(n, freq) {
	nind <- sum(n)
	inds <- 1:nind
	y <- rep(freq, times=n)
	z <- rep(inds,y)
	z
}

#' @export
genToFreq <- function(z) {
	x <- table(table(z))
	m <- data.frame(m=1:max(names(x)))
	x <- data.frame(x); colnames(x) <- c("m","freq")
	f <- merge(m, x, by.y=1, all=TRUE)
	f[is.na(f)] <- 0
	f$freq
}
