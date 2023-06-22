# Chao's 1987 estimator (based on mark-recapture, but commonly applied to species accumulation)

#' @export
chao.est <- function(capt.freq) {

	if (length(capt.freq) < 2) stop("At least 2 capture frquencies needed")

	sobs <- sum(capt.freq)
	a <- capt.freq[1]
	b <- capt.freq[2]

	if (sum(a,b) > sobs) stop("a + b can not exceed total number of observations")

	stot <- sobs + ((a^2) / (2*b))
	var.stot <- b*( 0.25*((a/b)^4) + ((a/b)^3) + 0.5*((a/b)^2) )
	sd.stot <- sqrt(var.stot)
	ll <- stot - 1.96*sqrt(var.stot)
	ul <- stot + 1.96*sqrt(var.stot)
	c(n=stot, var=var.stot, sd=sd.stot, ll=ll, ul=ul)
}
