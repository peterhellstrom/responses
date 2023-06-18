################################################################################
# POPULATION SIZE ESTIMATION BY USING ACCUMULATION CURVES
# Code for R
# By Peter Hellström, Department of Zoology, Stockholm University
# e-mail: peter.hellstrom@zoologi.su.se
# Last updated: September 1 2011
################################################################################

# Suggested updates:
# Sample from different distributions: assign probs. for Poisson & negative binomial. Check package sampling
# Enter recapture probability, and sample according to that.
# Parametric bootstrap sampling?
# Add output of variation for simulations

# IMPORTANT: Simulation part is not working at the moment, code is under re-development.
# Coding, try to remove excessive for-loops in simulation part, and clear up the code.
# Write a specific function that generates data, one that fits the models,
# and a couple of summary and plot-functions.
# How do for-loops compare with sapply and lapply? Not much difference on huge datasets!

################################################################################
# Required packages:
require(bbmle)
require(MASS)
require(nlme)
################################################################################
# This R-code uses three estimation functions that have been used in the literature

Kohn <- function(x,a,b) a*x/(b+x) # Michaelis-Menten
Eggert <- function(x,a,b) a*(1-exp(b*x)) # Two-parameter exponential
Chessel <- function(x,a) a-(a*(1-(1/a))^x)

Kohn.gr <- deriv(~ a*x / (b + x), c("a","b"), function(a,b,x) {})
Eggert.gr <- deriv(~ a*(1-exp(b*x)), c("a","b"), function(a,b,x) {})
Chessel.gr <- deriv(~ a-(a*(1-(1/a))^x), c("a"), function(a,x) {})

# Note that both the Michaelis-Menten equation and the two-parameter exponential
# equations are available in R as SSmicmen resp. SSasympOrig.
################################################################################
# REFERENCES:
# Kohn, M.H., York, E.C., Kamradt, D.A., Haught, G., Sauvajot, R.M., & Wayne, R.K. (1999) Estimating population size by genotyping faeces. Proceedings of the Royal Society B: Biological Sciences, 266, 657-663.
# Eggert, L.S., Eggert, J.A., & Woodruff, D.S. (2003) Estimating population sizes for elusive animals: the forest elephants of Kakum National Park, Ghana. Molecular Ecology, 12, 1389-1402.
# Frantz, A.C., Schaul, M., Pope, L.C., Fack, F., Schley, L., Muller, C.P., & Roper, T.J. (2004) Estimating population size by genotyping remotely plucked hair: the Eurasian badger. Journal of Applied Ecology, 41, 985-995.
# Bellemain, E., Swenson, J.E., Tallmon, D., Brunberg, S., & Taberlet, P. (2005) Estimating population size of elusive animals with DNA from hunter-collected faeces: four methods of brown bears. Conservation Biology, 19, 150-161.
# Petit, E. & Valiere, N. (2006) Estimating population size with noninvasive capture-mark-recapture data. Conservation Biology, 20, 1062-1073.
# Frantz, A.C. & Roper, T.J. (2006) Simulations to assess the performance of different rarefaction methods in estimating population size using small datasets. Conservation Genetics, 7, 315-318.
################################################################################

################################################################################
# Control settings for nls
nlc <- nls.control(maxiter = 10000, tol=0.000001, minFactor=0.000001)

################################################################################
# Calculate the cumulative number of unique value in a vector
cumunique <- function(x) cumsum(ifelse(!duplicated(x),1,0))

################################################################################
# Function for generating starting parameters for nls.
# Input x = samples and y = genotypes, unless iniMethod is grouped.
# Defaults to iniMethod="nls", available options are
# iniMethod = c("nls", "rlm", "lm")
# nls gives best values, then followed by rlm and lm.
# However, computation time is almost twice as long for bls and rlm compared with lm.
# It is also possible to use manually entered starting values in vector format.

get.ini.mm <- function(x, y=NULL, iniMethod="nls") {
	
	if (is.numeric(iniMethod) == TRUE) out <- iniMethod # If input is vector with values
	
	if (is.numeric(iniMethod) == FALSE) {
		
		if (iniMethod == "grouped") {
			# grouped can only be used if the function resampl.data() has been used to create an object
			# with resampled data sets
			yv.mean <- tapply(x[,"genotypes"], x[,"samples"], mean)
			data.set <- list(x = unique(x[,"samples"]), y = yv.mean)
			Kohn.ini <- try(getInitial(y ~ SSmicmen(x, a, b), data=data.set, control=nlc), silent=TRUE)
			out <- as.vector(Kohn.ini)
		}
		if (iniMethod == "nls") {
			# Use SSmicmen to get starting values
			# Create input data for nls
			# Remove first observation, for the case if y=log(y). Otherwise code crashes when y=0.
			inds <- which(y==0)
			if (length(inds) != 0) data.set <- list(x=x[-inds], y=y[-inds])
			if (length(inds) == 0) data.set <- list(x=x, y=y)
			# Use the built-in function SSmicmen to get initial parameter estimates for the "Kohn-model"
			Kohn.ini <- try(getInitial(y ~ SSmicmen(x, a, b), data=data.set, control=nlc), silent=TRUE)
			out <- as.vector(Kohn.ini)
		}
		if (iniMethod == "rlm") {
			# Use Lineweaver-Burk method and robust regression
			fm.lwb <- rlm(I(1/y) ~ I(1/x))
			coefs <- coef(fm.lwb)
			a <- 1/coefs[1]
			b <- a*coefs[2]
			out <- as.numeric(c(a,b))
		}
		if (iniMethod == "lm") {
			# Use Lineweaver-Burk method, ordinary least squares
			fm.lwb <- lm(I(1/y) ~ I(1/x))
			coefs <- coef(fm.lwb)
			a <- 1/coefs[1]
			b <- a*coefs[2]
			out <- as.numeric(c(a,b))
		}
	}
	out
}
################################################################################
# Input x = samples and y = genotypes, nS = length(samples), IC = Information criteria (returns a table with AICc and Akaike weights if TRUE),
# iniMethod = method for generation of start values
acc.fit <- function(x, y, nS, IC=TRUE, iniMethod="nls") {
	
	# Use get.ini.mm() to get starting values
	inits <- get.ini.mm(x=x, y=y, iniMethod=iniMethod)
	
	if (length(inits) != 2) message("Initial parameter estimation failed!", "\n")
	a <- inits[1]; b <- inits[2]

	# Inital estimates for Eggert and Chessel obtained by trial & error, will hopefully work for most datasets.
	# This section could be improved!
	mKohn <- try(nls(y ~ Kohn.gr(a,b,x), start=list(a=a, b=b), control=nlc), silent=TRUE)
	mEggert <- try(nls(y ~ Eggert.gr(a,b,x), start=list(a=a/2, b=-b/1000), control=nlc), silent=TRUE)
	mChessel <- try(nls(y ~ Chessel.gr(a,x), start=list(a=a/2), control=nlc), silent=TRUE)

	# Make a table that summarizes the parameter estimates for the three models
	A <- rbind(coef(mKohn), coef(mEggert), c(coef(mChessel),NA))
	rownames(A) <- c("mKohn","mEggert","mChessel")
		
	# Model selection
	# Returns a table with model selection parameters (uses the function ICtab from package bbmle)
	if (IC == TRUE) {
		ms.tab <- ICtab(mKohn,mEggert,mChessel, type="AICc", nobs=nS, weights=TRUE, delta=TRUE, base=TRUE, sort=FALSE)
	}
	if (IC == FALSE) {
		ms.tab <- NULL
	}
		
	list(A=A, ms.tab=ms.tab)
	}

################################################################################
# This function fits all three models and returns a table with the parameter estimates.
acc.curve <- function(samples, genotypes, shuffle=TRUE, iniMethod="nls", Log=FALSE) {
	
	nS <- length(samples)
	nUG <- length(unique(genotypes))
	
	x <- 1:nS
	z <- rep(1:nUG, as.vector(table(genotypes)))
	if (shuffle == TRUE) {
		# Shuffle the dataset once before analysis
		z <- sample(z, nS, replace=FALSE, prob=NULL)
		}
	# Calculate the cumulative number of unique genotypes in the sample z
	y <- cumunique(z)
	if (Log==TRUE) y <- log(y)

	fitm <- acc.fit(x, y, nS=nS, iniMethod=iniMethod)
	
		out <- list(
			"SampleSize" = nS,
			"NoUniqueGenotypes" = nUG,
			"ParameterEstimates" = fitm$A,
			"msel" = fitm$ms.tab,
			"Samples" = x,
			"Genotypes" = y)
			
		# Print output
		out
}

################################################################################
plot.acc.curve <- function(obj, plot.type="normal", xlim=NULL, addPoints=TRUE) {

	x <- obj$Samples
	y <- obj$Genotypes
	A <- obj$ParameterEstimates
	
	if (is.null(xlim)) {
		# Make a plot showing a) the accumulation curve and b) the three different fitted models
		# Extended plot type means
		if (plot.type=="extended") { xlims <- c(0,A[1,1]*A[1,2]-A[1,2]); ylims <- c(0, A[1,1]) }
		else if (plot.type=="normal") { xlims <- c(0,max(x)); ylims <- c(0,max(y)) }
		else stop ("The plot type you specified is not valid!") 
		}
	
	if (!is.null(xlim)) {
		# Make a plot showing a) the accumulation curve and b) the three different fitted models
		# Extended plot type means
		if (plot.type=="extended") { xlims <- xlim; ylims <- c(0, A[1,1]) }
		else if (plot.type=="normal") { xlims <- xlim; ylims <- c(0,max(y)) }
		else stop ("The plot type you specified is not valid!") 
		}
		
		dev.new(width=12,height=6)
		par(mfrow=c(1,2))

		plot(c(0,x), c(0,y), pch=16, bty="l", xlab="Sample number",
			ylab="Unique genotypes", type="l", lwd=2, lty=1,
			ylim=ylims, xlim=c(0,max(x)),
			main="Accumulation curve", font.lab=2, las=1)
			points(c(0,x), c(0,y), col="red", pch=16, cex=0.7)

		plot(x=c(0,xlims), c(0,ylims), pch=16, bty="l", xlab="Sample number",
			ylab="Unique genotypes", type="n", lwd=2, lty=1,
			ylim=ylims, xlim=xlims,
			main="Model estimates", font.lab=2, las=1)

			if (addPoints == TRUE) points(c(0,x), c(0,y), col="red", pch=16, cex=0.7)
			curve(A[1,1]*x/(A[1,2]+x), from=0, to=xlims[2], add=T, type="l", lty=1, lwd=2, col="black")
			curve(A[2,1]*(1-exp(A[2,2]*x)), add=T, type="l", lty=2, lwd=2, col="green")
			curve(A[3,1]-(A[3,1]*(1-(1/A[3,1]))^x), add=T, type="l", lty=3, lwd=2, col="blue")
			abline(v=length(x), lty=2) # Draw vertical line at number of samples

			legend("topleft",c("Kohn","Eggert","Chessel"), col=c("black","green","blue"), lwd=c(2,2,2), lty=c(1,2,3), bty="n", cex=0.8)
		
			par(mfrow=c(1,1))
}

################################################################################
# This function randomizes the input data and runs the analysis a specified number of times.
# It uses the acc.curve() function

resampl.data <- function(samples, genotypes, resamplings, replace=FALSE, Log=FALSE) {
	
	nS <- length(samples)
	nUG <- length(unique(genotypes))
	nR <- resamplings
	x <- 1:nS
	# Create a matrix with datasets
	z <- sapply(1:resamplings, function(i) sample(genotypes, nS, replace=replace, prob=NULL) )
	# Calculate the cumulative number of unique genotypes in the sample z
	y <- sapply(1:resamplings, function(i) cumunique(z[,i]))
	if (Log == TRUE) y <- log(y)
	rownames(y) <- x
	colnames(y) <- paste("Smp", 1:resamplings, sep="")
	
	yLong <- data.frame(
				genotypes = c(y),
				samples = rep(x, times=resamplings),
				resample = as.factor(rep(colnames(y), each=nS)))
	yGr <- groupedData(genotypes ~ samples|resample, data=yLong)

	out <- list(
		nS = nS,
		nUG = nUG,
		resamplings = nR,
		Log = Log,
		x = x,
		z = z,
		y = y,
		yGr = yGr)
	out
}
################################################################################
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
################################################################################
# Fit nonlinear mixed models
acc.curve.mixed <- function(input.data) {
	
	nS <- input.data$nS
	nUG <- input.data$nUG
	rS <- input.data$resamplings
	x <- input.data$x
	yGr <- input.data$yGr
	
	coefs <- get.ini.mm(x=input.data$yGr, iniMethod="grouped")
	a <-coefs[1]; b <-coefs[2]

	# Fit nonlinear mixed models:
	mKohn <- nlme(genotypes ~ Kohn.gr(a=a, b=b, x=samples),
		fixed = a + b ~ 1, random = a + b ~ 1, start = c(a=a, b=b), data = yGr)

	mEggert <- nlme(genotypes ~ Eggert.gr(a=a, b=b, x=samples),
		fixed = a + b ~ 1, random = a + b ~ 1, start = c(a=a/2, b= -b/1000), data = yGr)

	mChessel <- nlme(genotypes ~ Chessel.gr(a=a, x=samples),
		fixed = a ~ 1, random = a ~ 1, start = c(a=a/2), data = yGr)
	
	# Extract data from fits
	regr <- cbind(coef(mKohn), coef(mEggert), coef(mChessel))
	regr <- regr[,c(1,3,5,2,4)]
	colnames(regr) <- c("a.Kohn", "a.Eggert", "a.Chessel", "b.Kohn", "b.Eggert")
			
	# Model selection
	# Returns a table with model selection parameters (uses the function ICtab from package bbmle)
	ms.tab <- ICtab(mKohn, mEggert, mChessel, type="AICc", nobs=nS, weights=TRUE, delta=TRUE, base=TRUE, sort=FALSE)
	
	if (input.data$Log == TRUE) regr <- exp(regr)
	
	# Extract summary statistics from the resampling parameter estimates
	A <- apply(regr, 2, summary)
	sds <- apply(regr, 2, sd)
	qts <- apply(regr, 2, quantile,c(0.025,0.975))
	skews <- apply(regr, 2, skew); kurts <- apply(regr, 2, kurtosis)
	A <- rbind(A, SD=sds, qts, Skew=skews, Kurtosis=kurts)
	
	out <- list(
		"SampleSize" = nS,
		"UniqueGenotypes" = nUG,
		"Resamplings" = rS,
		"ParameterSummary" = A,
		"ModelSelection" = 	ms.tab,
		"Coefficients" = regr,
		"ranef" = list("Kohn" = ranef(mKohn), "Eggert" = ranef(mEggert), "Chessel" = ranef(mChessel)),
		"fixef" = list("Kohn" = fixef(mKohn), "Eggert" = fixef(mEggert), "Chessel" = fixef(mChessel)),
		"summary" = list("Kohn" = summary(mKohn), "Eggert" = summary(mEggert), "Chessel" = summary(mChessel)),
		"models" = list("Kohn" = mKohn, "Eggert" = mEggert, "Chessel" = mChessel)
		)
		
	out
}
################################################################################
# Plot the population size estimates (min,max & mean) from acc.curve.resamp(), 
# and also the distribution of the population size estimate 'a'.
# Requires a stored object fitted with acc.curve.resamp().

plot.curve <- function(obj, model, ylim=NULL) {

	ests <- obj$Coefficients[,1:5]
	
	ests.mean <- apply(ests, 2, mean)
	ests.median <- apply(ests, 2, median)
	
	# ests.min <- apply(ests, 2, min)
	# ests.max <- apply(ests, 2, max)
	ests.min <- sapply(1:ncol(ests), function(i) quantile(ests[,i],c(0.025)))
	ests.max <- sapply(1:ncol(ests), function(i) quantile(ests[,i],c(0.975)))
	
	if (model == "Kohn") inds <- c(1,4)
	if (model == "Eggert") inds <- c(2,5)
	if (model == "Chessel") inds <- 3

	xlims <- c(0, obj$SampleSize)
	if (is.null(ylim)) ylims <- c(0, ests.max[inds][1])
	if(!is.null(ylim)) ylims <- ylim
	
	dev.new(width=7,height=6)
	plot(x=xlims, y=ylims, type="n", xlab="Sample number", ylab="Discovered genotypes", 
		font.lab=2, las=1, bty="l", main=paste("Population size based on ", model, " equation", sep=""))
	
	abline(h=ests.min[inds][1], lty=2, col=2)
	abline(h=ests.max[inds][1], lty=2, col=2)
	abline(h=ests.median[inds][1], lty=1, col=2)
	
	if (model != "Chessel") {
		str.min <- parse(file="", text=paste("curve(",model, "(x=x, a=ests.min[inds][1], b=ests.min[inds][2]), n=101, lty=2, add=T)", sep=""))
		str.max <- parse(file="", text=paste("curve(",model, "(x=x, a=ests.max[inds][1], b=ests.max[inds][2]), n=101, lty=2, add=T)", sep=""))
		str.median <- parse(file="", text=paste("curve(",model, "(x=x, a=ests.median[inds][1], b=ests.median[inds][2]), n=101, lty=1, add=T)", sep=""))
		
		eval(str.min)
		eval(str.max)
		eval(str.median)
		}
		
	if (model == "Chessel") {
		curve(Chessel(x, a=ests.min[inds]), from=0, to=obj$SampleSize, n=101, lty=2, add=T)
		curve(Chessel(x, a=ests.max[inds]), from=0, to=obj$SampleSize, n=101, lty=2, add=T)
		curve(Chessel(x, a=ests.median[inds]), from=0, to=obj$SampleSize, n=101, lty=1, add=T)
	}

}

################################################################################
plot.hist <- function(obj, model, breaks=30, plot.type="histogram", drawN=TRUE) {

	rs <- obj$Resamplings
	nS <- obj$SampleSize
	ests <- obj$Coefficients[,1:3]
	colnames(ests) <- c("Kohn","Eggert","Chessel")

	dev.new(width=8,height=6)
	# Plot distribution of estimates
	
	if (plot.type == "histogram") {
	# Histogram
		store.hist <- hist(ests[,model], freq=F, font.lab=2, las=1,
					breaks=breaks, col="lightgray", xlab="Estimate", ylab="Density", 
					main=paste("Population size estimate,", model))
		if (drawN == TRUE) {
			# Overlay normal distribution
			.x <- seq(min(store.hist$breaks), max(store.hist$breaks), 0.1)
			lines(.x, dnorm(.x, mean=mean(ests[,model]), sd=sd(ests[,model])), type="l", lwd=2, col=4)
			}
		}
	if (plot.type == "density") {
		plot(density(ests[,model]), lwd=2, col=4, las=1, font.lab=2, main=paste("Population size estimate,", model))
		}
		
	# Vertical lines for max, min & median estimates
	abline(v=quantile(ests[,model],0.025), lty=2, lwd=2, col="red")
	abline(v=quantile(ests[,model],0.975), lty=2, lwd=2, col="red")
	abline(v=median(ests[,model]), lty=1, lwd=2, col="black")
	
	# Text showing mean value
	title(sub=paste("Median =", round(median(ests[,model]),2), ",", rs, "resamplings"), cex=0.8, font=3)
	intrvl <- round(as.numeric(quantile(ests[,model], c(0.025, 0.975))),2)
	mtext(text = paste("95% Quantile interval =", intrvl[1],"-",intrvl[2]), side=3, line=0.5)
}

################################################################################
# Make a plot that compares means and q-Intervals for the three different models
# Requires a stored object fitted with acc.curve.resamp().

plot.est.comp <- function(obj) {

	A <- obj$ParameterSummary

	dev.new(width=8,height=6)
	# Set plotting region and plot mean values
	plot(A["Mean",1:3], xlim=c(0.75, 3.25), ylim=c(min(A["Min.",1:3]),max(A["97.5%",1:3])),
		pch=15, cex=1.75, col="red",
		ylab="Population size", xlab="Estimation function", 
		bty="l", xaxt="n", las=1, font.lab=2, font=2, main="Means, medians & 95% quantile intervals")
	# Add median values to the plot
		points(x=c(1,2,3), y=A["Median",1:3], pch=16, cex=1.75, col="blue")
	# Draw the x-axis
		axis(1, c(1,2,3), tcl=0.3, labels=c("Kohn","Eggert","Chessel"), font=2)
	# Draw the 95% quantile interval
		arrows(c(1,2,3), A["2.5%",1:3], c(1,2,3), A["97.5%",1:3], length=0.11, angle=90, code=3, lwd=2)
	# Add legend
	legend("topright", c("Mean", "Median"), col=c(2,4), pch=c(15,16), cex=1, bty="n")
	}

################################################################################
# Plot correlation between parameter estimates for Kohn and Eggert models:
plot.par.cor <- function(obj, ...) {
	dev.new(width=12,height=6)
	par(mfrow=c(1,2))
	plot(obj$Coefficients[,"b.Kohn"], obj$Coefficients[,"a.Kohn"], cex=0.7, 
	xlab="Half-saturation constant", ylab="Asymptote", font.lab=2, las=1, main="Kohn")
	abline(h=obj$ParameterSummary["Mean","a.Kohn"], col=2, lty=2, ...)

	plot(obj$Coefficients[,"b.Eggert"], obj$Coefficients[,"a.Eggert"], cex=0.7,
	xlab="Natural logarithm of the rate constant", ylab="Asymptote", font.lab=2, las=1, main="Eggert")
	abline(h=obj$ParameterSummary["Mean","a.Eggert"], col=2, lty=2, ...)
	par(mfrow=c(1,1))
}
	
################################################################################
# This function computes the estimate of skewness of the distribution
# of x, following Box 6.2 in Sokal and Rohlf
skew <- function(x) {
	m <- mean(x)
	s <- sd(x)
	y <- x - m
	sumy3 <- sum(y^3)
	n <- length(x)
	( n/((n-1)*(n-2)) ) * sumy3/s^3
	}
################################################################################
# This function computes the estimate of kurtosis of the distribution
# of x, following Box 6.2 in Sokal and Rohlf
kurtosis <- function(x) {
	m <- mean(x)
	s <- sd(x)
	y <- x - m
	sumy4 <- sum(y^4)
	n <- length(x)
	( (n+1)*n/((n-1)*(n-2)*(n-3)) ) * sumy4/s^4 - 3*(n-1)^2/((n-2)*(n-3))
}
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
################################################################################
# PART 2:
# Simulation functions
################################################################################
################################################################################

# Restructure code like this instead.
# sim.gen.data <- function(genotypes, samples, resamplings) ()
# sim.run.data <- function(input.data) ()
# sim.summary <- function(run.data) ()
# sim.plot <- function(sim.summary) ()

acc.curve.sim <- function(genotypes, samples, resamplings) {
	
	# Inputs: genotypes (list with true population sizes), samples, resamplings

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )

	# Vectors to hold output (g is for global)
	avg.output <- vector("list",length=length(samples))
	full.output <- vector("list",length=length(samples))
	g.avg.output <- vector("list",length(genotypes))
	g.full.output <- vector("list",length(genotypes))
	avg.prop.found <- matrix(NA,ncol=length(genotypes),nrow=length(samples))
	colnames(avg.prop.found) <- pop.vec
	rownames(avg.prop.found) <- samples 
	full.prop.found <- matrix(NA,ncol=resamplings,nrow=length(samples))
	g.full.prop.found <- vector("list",length(genotypes))   

	# Outer loop - loops over length of true pop size (input genotypes)
	# Inner loop - sample size

	for (runs in 1:length(genotypes)) { 

		for (res in 1:length(samples)) {

		# Create data sets
		x <- 1:samples[res]
		y <- matrix(NA,ncol=resamplings,nrow=samples[res])
			for (j in 1:resamplings) y[,j]<- sample(genotypes[[runs]],samples[res],replace=TRUE)

		z <- matrix(NA,ncol=resamplings,nrow=samples[res]) 
			for (k in 1:resamplings) z[,k] <- rep(1:length(unique(y[,k])),as.vector(table(y[,k])) )
			z.mat <- matrix(rep(z,length(resamplings)),nrow=samples[res],ncol=resamplings,byrow=FALSE)
			shuffle <- apply(z.mat,2,sample,replace=FALSE)

		data.set <- matrix(NA,nrow=dim(shuffle)[1],ncol=dim(shuffle)[2])
			for (l in 1:ncol(shuffle)) {
				for (m in 1:nrow(shuffle)) data.set[m,l] <- length(unique(shuffle[1:m,l]) ) }

		
		# Extract proportion of found genotypes (averaged over resamplings)
		a <- max(genotypes[[runs]])
		b <- nrow(data.set)
		avg.prop.found[res,runs] <- ("Mean prop found"=as.vector(summary(data.set[b,])[4])/a )

		# Extract proportion found for each individual run
		full.prop.found[res,] <- data.set[b,]/a

		# Non-linear curve fitting
		kohn <- lapply(1:resamplings, function(i) try(summary(nls(y~SSmicmen(x,a,b),data=list(x=x,y=data.set[,i]))) ) )

		kohn <- kohn[sapply(kohn, function(x) !inherits(x, "try-error"))]
		kohn.l <- length(kohn)
		kohn.est <- t(sapply(1:kohn.l, function(i) kohn[[i]]$parameters[1,]))

		eggert <- lapply(1:resamplings, function(i) try(summary(nls(y~SSasympOrig(x,a,b),data=list(x=x,y=data.set[,i]))) ) )
		eggert <- eggert[sapply(eggert, function(i) !inherits(i, "try-error"))]
		eggert.l <- length(eggert)
		eggert.est <- t(sapply(1:eggert.l, function(i) eggert[[i]]$coefficients[1,]))

		chessel <- lapply(1:resamplings, function(i) try(summary(nls(y~Chessel(x,a),start=list(a=median(eggert.est[,1])),data=list(x=x,y=data.set[,i])) ) ) )
		chessel <- chessel[sapply(chessel, function(i) !inherits(i, "try-error"))]
		chessel.l <- length(chessel)
		chessel.est <- t(sapply(1:chessel.l, function(i) chessel[[i]]$coefficients[1,]))

	# Create output with average values (stored in avg.output)
	quants <- function (x) as.vector(quantile(x,c(0.025,0.975) ) ) # Get 95% quantile interval
	# Gather quantile intervals in a matrix
	qi <- matrix(c(quants(kohn.est[,1]),quants(eggert.est[,1]),quants(chessel.est[,1])),ncol=3,nrow=2)
	# Gather number of succesful runs
	coef.l <- as.vector(c(length(kohn.est[,1]),length(eggert.est[,1]),length(chessel.est[,1])))
	# Combine parameter estimates & intervals in one output
	avg.out <- rbind(cbind(as.vector(summary(kohn.est[,1])), as.vector(summary(eggert.est[,1])), as.vector(summary(chessel.est[,1])) ),qi,coef.l)
	# Set row- & colnames
	rownames(avg.out) <- c("Min","1st Qu.","Median","Mean","3rd Qu.","Max.","2.5% Qu.","97.5% Qu.","N")
	colnames(avg.out) <- c("kohn","eggert","chessel")
	# Write avg.output to list-object
	avg.output[[res]] <- avg.out
	full.output[[res]] <- list(kohn.est=kohn.est,eggert.est=eggert.est,chessel.est=chessel.est)
	}
	g.avg.output[[runs]] <- avg.output
	g.full.output[[runs]] <- full.output
	g.full.prop.found[[runs]] <- full.prop.found
	}
	list(average=g.avg.output,full=g.full.output,AvgPropFound=avg.prop.found,FullPropFound=g.full.prop.found)
}

################################################################################

# Data extraction and summarization of simulated datasets:
# After the looping-procedure: extract the averaged data generated by the output
# from the acc.curve.sim-function.
# Requires a stored object with all data generated by the acc.curve.sim-function

acc.curve.xtr.avg <- function(genotypes, samples, object, statistic) {

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )

	# Extract the population size estimates
	avg.est <- vector("list",length=length(genotypes)) 
	# 4 for mean, 3 for median
	if (statistic=="mean")
		stat <- 4
	else if (statistic=="median")
		stat <- 3
	else stop ("Not a valid statistic parameter!")

	for (j in 1:length(genotypes)) {
		avg.dat <- sapply(1:length(samples), function(i) as.vector(object[[1]][[j]][[i]][stat,]) )
		avg.est[[j]] <- t(avg.dat)
		}

	avg <- function(genotypes, model) {
	dat <- sapply(1:length(genotypes), function (i) avg.est[[i]][,model])
	colnames(dat) <- pop.vec
	rownames(dat) <- samples
	dat
	}

	kohn.avg <- avg(genotypes,1) 
	eggert.avg <- avg(genotypes,2)
	chessel.avg <- avg(genotypes,3)

	# Extract quantile intervals
	qi.est <- vector("list",length=length(genotypes))

	for (j in 1:length(genotypes)) {
		qi.dat <- sapply(1:length(samples), function(i) as.vector(object[[1]][[j]][[i]][c(7,8),]) )
		qi.est[[j]] <- t(qi.dat)
		}

	qis <- function(genotypes,model) {
		if (model==1) ind <- c(1,2)
		if (model==2) ind <- c(3,4)
		if (model==3) ind <- c(5,6)

	dat <- lapply(1:length(genotypes), function (i) qi.est[[i]][,ind])
		for (i in 1:length(genotypes)) colnames(dat[[i]]) <- c("2.5%","97.5%")
		for (i in 1:length(genotypes)) rownames(dat[[i]]) <- samples
	dat
	}

	kohn.qi <- qis(genotypes,1)
	eggert.qi <- qis(genotypes,2)
	chessel.qi <- qis(genotypes,3)

	# Estimate bias, each function at the time

	TruePop <- matrix(rep(pop.vec,length(samples)),nrow=length(samples),ncol=length(genotypes),byrow=T )

	# Bias function
	# model: 1 = Kohn, 2 = Eggert, 3 = Chessel
	bias.est <- function(input, TruePop) {
	bias.mat <- (input-TruePop) / TruePop
	dat <- 100*bias.mat
	colnames(dat) <- pop.vec
	rownames(dat) <- samples
	dat
	}

	kohn.bias <- bias.est(kohn.avg,TruePop)
	eggert.bias <- bias.est(eggert.avg,TruePop)
	chessel.bias <- bias.est(chessel.avg,TruePop)

	# Create output
	list(
		kohn = list(samples = samples, avg.est = kohn.avg, qi = kohn.qi, bias = kohn.bias, prop.found = object[[3]]),
		eggert = list(samples = samples, avg.est = eggert.avg, qi = eggert.qi, bias = eggert.bias, prop.found = object[[3]]),
		chessel = list(samples = samples, avg.est = chessel.avg, qi = chessel.qi, bias = chessel.bias, prop.found = object[[3]])
	)
}

################################################################################
# Plot of average values
# (Estimated pop size against sample size)

acc.sim.avg.plot <- function(genotypes, samples, object, model, type) {

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )

	model.lab <- model

	if (model=="kohn") model <- 1
	else if (model=="eggert") model <- 2
	else if (model=="chessel") model <- 3
	else stop ("Wrong model chosen!")

	if (type=="qi") {
		xlims <- max(object[[model]]$samples) 
		ylims <- max(unlist(object[[model]]$qi))

	dev.new(width=6,height=5)
	par(mar=c(5,4,3,1))
	plot(c(0,xlims),c(0,ylims),type="n",xlab="Sample size",ylab="Estimated population size",font.lab=2,las=1,main=model.lab)

	# Plot point estimates
	xv <- rep(object[[model]]$samples,times=length(genotypes) )
	yv <- as.vector(object[[model]]$avg.est)
	cols <- rep(1:length(genotypes), each=length(samples))
	points(xv,yv,pch=16,col=cols,cex=1.3)

	# Plot quantile intervals
	for (i in 1:length(genotypes)) {
		xvec <- samples
		low.arr <- object[[model]]$qi[[i]][,1]
		up.arr <- object[[model]]$qi[[i]][,2]
		arrows(xvec,low.arr,xvec,up.arr,angle=90,code=3,length=0.05) }

	# Plot true population size as horizontal line
	pop.vec <- colnames(object[[model]]$avg.est)
	for (i in 1:length(genotypes)) {
		abline(h=pop.vec[i],lty=2,col=c(0+i))
	}} 

	if (type=="point") {
		xlims <- max(object[[model]]$samples) 
		ylims <- max(unlist(object[[model]]$avg.est))

		dev.new(width=6,height=5)
		par(mar=c(5,4,3,1))
		plot(c(0,xlims),c(0,ylims),type="n",xlab="Sample size",ylab="Estimated population size",font.lab=2,las=1,main=model.lab)

		# Plot point estimates
		xv <- rep(object[[model]]$samples,times=length(genotypes) )
		yv <- as.vector(object[[model]]$avg.est)
		cols <- rep(1:length(genotypes), each=length(samples))
		points(xv,yv,pch=16,col=cols,cex=1.3)

		# Plot true population size as horizontal line
		pop.vec <- colnames(object[[model]]$avg.est)
		for (i in 1:length(genotypes)) {
		abline(h=pop.vec[i],lty=2,col=c(0+i))
		}}

	# Add legend
	legend("topright",legend=pop.vec, col=c(1:length(genotypes)),pch=16,cex=1.3 ) 
}


################################################################################
# Plot of average values
# Estimated relative bias against sample size)

acc.sim.bias.plot <- function(genotypes, samples, object, model) {

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )

	model.lab <- model

	if (model=="kohn") model <- 1
	else if (model=="eggert") model <- 2
	else if (model=="chessel") model <- 3
	else stop ("Wrong model chosen!")

	xlims <- max(object[[model]]$samples)
	ylims.lo <- min(unlist(object[[model]]$bias))
	ylims.up <- max(unlist(object[[model]]$bias))

	dev.new(width=6,height=5)
	par(mar=c(5,4,3,1))
	plot(c(0,xlims),c(ylims.lo,ylims.up),type="n",xlab="Sample size",ylab="Relative bias",font.lab=2,las=1,main=model.lab)

	# Plot point estimates
	xv <- rep(object[[model]]$samples,times=length(genotypes) )
	yv <- as.vector(object[[model]]$bias)
	cols <- rep(1:length(genotypes), each=length(samples))
	points(xv,yv,pch=16,col=cols,cex=1.3)

	# Add legend
	legend("topright",legend=pop.vec, col=c(1:length(genotypes)),pch=16,cex=1.3 )
	}
################################################################################

acc.sim.bias.ratio.plot <- function(genotypes, samples, object, model) {

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )
	.x <- sapply(1:length(genotypes), function(i) samples/pop.vec[i])
 
	model.lab <- model

	if (model=="kohn") model <- 1
	else if (model=="eggert") model <- 2
	else if (model=="chessel") model <- 3
	else stop ("Wrong model chosen!")

	xlims <- max(.x)
	ylims.lo <- min(unlist(object[[model]]$bias))
	ylims.up <- max(unlist(object[[model]]$bias))

	dev.new(width=6,height=5)
	par(mar=c(5,4,3,1))
	plot(c(0,xlims),c(ylims.lo,ylims.up),type="n",xlab="Sample size / True population size",ylab="Relative bias",font.lab=2,las=1,main=model.lab)

	# Plot point estimates
	xv <- as.vector(.x)
	yv <- as.vector(object[[model]]$bias)
	cols <- rep(1:length(genotypes), each=length(samples))
	points(xv,yv,pch=16,col=cols,cex=1.3)

	# Add legend
	legend("topright",legend=pop.vec, col=c(1:length(genotypes)),pch=16,cex=1.3 )
}

################################################################################
################################################################################

################################################################################
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
################################################################################
