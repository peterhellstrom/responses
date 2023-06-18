n <- 200
nrepl <- 100

prob <- runif(200)
prob

chao.table3 <- function(n, t, prob, nrepl) {
	
	x <- lapply(1:nrepl, function(i) {
		x <- replicate(sample(1:n, prob=prob, replace=TRUE), n=t)
		table(x)
		})
	
	# estimates
	estimate <-
	se <- 
	coverage <-
	
}


y=sample(1:3,15,1)

prob=matrix(runif(45),15) 

prob[cbind(1:length(y), y)]



################################################################################
# Chao's 1987 estimator (based on mark-recapture, but commonly applied to species accumulation)

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

genToFreq <- function(z) {
	x <- table(table(z))
	m <- data.frame(m=1:max(names(x)))
	x <- data.frame(x); colnames(x) <- c("m","freq")
	f <- merge(m, x, by.y=1, all=TRUE)
	f[is.na(f)] <- 0
	f$freq
}

FreqTogen <- function(n, freq) {
	nind <- sum(n)
	inds <- 1:nind
	y <- rep(freq, times=n)
	z <- rep(inds,y)
	z
}
################################################################################

