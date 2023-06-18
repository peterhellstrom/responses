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
