ProjectSurvival <- function(N0, T, st) {
	# N0 = initial number of individuals
	# T = total number of transitions ("years"-1)
	# st = time specific survival (vector of same length as T)
	if (length(st) != T) stop("st must be a vector with length T")
	out = matrix(ncol=T+1, nrow=length(N0))
	out[,1] <- N0
	for (t in 1:T) {
		N.surv = rbinom(n=N0, size=N0, prob=st[t])
		out[,t+1] <- N.surv
		N0 = N.surv
	}
	colnames(out) <- paste("t",1:ncol(out),sep="")
	rownames(out) <- paste("b",1:nrow(out),sep="")
	out
}
