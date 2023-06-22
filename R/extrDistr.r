# Calculates IC and IC weights, option order=TRUE/FALSE.
# order=TRUE returns the table with candidate models, in order according to AICc (from smallest to largest).

#' @export
extrDistr <- function(x, order=FALSE) {
	.x <- extrDistrAIC(x)
	distname <- sapply(1:length(x), function(i) x[[i]]$distname)
	colnames(.x) <- distname
	rownames(.x) <- c("n","k","logLik","AIC","AICc","BIC")
	.x <- t(.x)
	#.x
	delta <- .x[,"AICc"] - min(.x[,"AICc"])
	L <- exp(-delta/2) # likelihoods of models
	w <-  L / sum(L) # Akaike weights
	.x <- data.frame(.x,delta,w)

	if (order==TRUE) {
		inds <- order(delta)
		.x <- .x[inds,]
		}
	.x
}


# Input is a list with objects fitted with fitdist from fitdistrplus
# The input list must be created manually.
# This function is an intermediate function, and calculates various information criteria
# It's intended use is ONLY from within extrDistr.
#' @export
extrDistrAIC <- function(x) {

	if (class(x)!="list") stop("Input must be a list")

	sapply(1:length(x), function(i) {

		n <- length(x[[i]]$data)
		k <- length(x[[i]]$estimate)

		loglik <- x[[i]]$loglik

		aic <- 2*k - 2*loglik
		# aic <- x[[i]]$aic
		aic.c <- -2*loglik + 2*k*(k+1)/(n-k-1)
		bic <- -2*loglik + k*log(n)

		.x <- c(n,k,loglik,aic,aic.c,bic)
	})
}
