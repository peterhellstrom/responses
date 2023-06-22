# Input x = samples and y = genotypes, nS = length(samples), IC = Information criteria (returns a table with AICc and Akaike weights if TRUE),
# iniMethod = method for generation of start values
#' @export
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

# This function fits all three models and returns a table with the parameter estimates.
#' @export
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
