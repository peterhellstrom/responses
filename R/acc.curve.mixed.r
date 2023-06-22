# Fit nonlinear mixed models
#' @export
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
