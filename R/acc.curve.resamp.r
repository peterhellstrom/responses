acc.curve.resamp <- function(input.data, iniMethod="nls", IC=TRUE, distmom=FALSE) {
	
	nS <- input.data$nS
	nUG <- input.data$nUG
	rS <- input.data$resamplings
	x <- input.data$x
	
	res <- lapply(1:rS, function(i) {
		try(acc.fit(x, input.data$y[,i], nS=nS, IC=IC, iniMethod=iniMethod), silent=TRUE) }) ####
	
	# Calculate number of succesful resamplings
	rows <- length(res[sapply(res, function(x) !inherits(x, "try-error"))])
	if (rows==0) message("No valid cases.", "\n")
	else
	
	# Create a list that holds output for succesful fits only.
	inds <- sapply(res, function(x) !inherits(x, "try-error"))
	fits <- res[inds]
	fail <- which(inds!=TRUE)
	
	# Extract data from succesful fits
	regr <- t(sapply(1:rows, function(i) unlist(fits[[i]]$A)))[,1:5]
	colnames(regr) <- c("a.Kohn", "a.Eggert", "a.Chessel", "b.Kohn", "b.Eggert")
	rownames(regr) <- which(inds==TRUE)
	
	# Create a vector that stores which model that provides the best fit as indicated by AICc.
	if (IC == TRUE) msel <- sapply(1:rows, function(i) which((fits[[i]]$ms.tab$dAICc)==0))
	if (IC == FALSE) msel <- NULL
	
	if (input.data$Log == TRUE) regr <- exp(regr)
	
	# Extract summary statistics from the resampling parameter estimates
	A <- apply(regr, 2, summary)
	sds <- apply(regr, 2, sd)
	qts <- apply(regr, 2, quantile,c(0.025,0.975))
	A <- rbind(A, SD=sds, qts)
	
	if (distmom == TRUE) {
		# Also calculate distribution parameters; skewness and kurtosis
		skews <- apply(regr, 2, skew); kurts <- apply(regr, 2, kurtosis)
		A <- rbind(A, Skew=skews, Kurtosis=kurts)
		}
	
	out <- list(
		"SampleSize" = nS,
		"UniqueGenotypes" = nUG, #### Why this one???
		"Resamplings" = rS,
		"SuccesfulResamplings" = rows,
		"FailedResamplings" = fail,
		"ParameterSummary" = A,
		"ModelSelection" = 	table(ifelse(msel==1, "mKohn", ifelse(msel==2, "mEggert", "mChessel"))),
		"Coefficients" = regr)
		
	out
}
