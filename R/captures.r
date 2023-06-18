# Sources:
# Otis DL, Burnham KP, White GC, Anderson DR (1978) Statistical inference from capture data on closed populations. Wildl Monogr 62:1-135 (p. 44-)
# White GC, Anderson DR, Burnham KP, Otis DL (1982) Capture-recapture and removal methods for sampling closed populations. Los Alamos National Laboratory, Los Alamos, New Mexico

captures <- function(N, sessions, c.prob) {

	E.n <- numeric(sessions)
	capt.n <- numeric(sessions)
	not.caught <- numeric(sessions)
	pop.start <- numeric(sessions)
	n.start <- N
	
	if (length(c.prob) == 1) c.prob <- rep(c.prob, sessions)
	if (length(c.prob) != sessions) stop("c.prob must be a single probability, or a vector of the same length as the number of sessions")

	for (j in 1:sessions) {
		E.n[j] <- N * ((1 - c.prob[j])^(j-1)) * c.prob[j] # This is probably incorrect if c.prob varies between sessions?
		capt.n[j] <- rbinom(1, n.start, c.prob[j])
		pop.start[j] <- n.start
		not.caught[j] <- n.start - capt.n[j]
		n.start <- not.caught[j]
	}

	mat <- matrix(ncol=5,nrow=sessions)
	mat[,1] <- 1:sessions
	mat[,2] <- pop.start
	mat[,3] <- E.n
	mat[,4] <- capt.n
	mat[,5] <- not.caught
	colnames(mat) <- c("Occasion", "Pop Size At Start", "Exp Captures", "Captures", "Not yet caught")
	mat
}
