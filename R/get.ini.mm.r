# Function for generating starting parameters for nls.
# Input x = samples and y = genotypes, unless iniMethod is grouped.
# Defaults to iniMethod="nls", available options are
# iniMethod = c("nls", "grouped", "rlm", "lm")
# nls gives best values, then followed by rlm and lm.
# However, computation time is almost twice as long for nls and rlm compared with lm.
# It is also possible to use manually entered starting values in vector format.

#' @export
get.ini.mm <- function(x, y=NULL, iniMethod=c("nls","grouped","rlm","lm")) {

	if (is.numeric(iniMethod) == TRUE) out <- iniMethod # If input is vector with values

	if (is.numeric(iniMethod) == FALSE) {

		iniMethod <- match.arg(iniMethod)

		if (iniMethod == "grouped") {
			# grouped can only be used if the function resampl.data() has been used to create an object
			# with resampled data sets
			yv.mean <- tapply(x[,"genotypes"], x[,"samples"], mean)
			data.set <- list(x = unique(x[,"samples"]), y = yv.mean)
			Kohn.ini <- try(getInitial(y ~ SSmicmen(x, a, b), data=data.set, control=nlc), silent=TRUE)
			out <- as.vector(Kohn.ini)
		}
		else if (iniMethod == "nls") {
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
		else if (iniMethod == "rlm") {
			# Use Lineweaver-Burk method and robust regression
			fm.lwb <- rlm(I(1/y) ~ I(1/x))
			coefs <- coef(fm.lwb)
			a <- 1/coefs[1]
			b <- a*coefs[2]
			out <- as.numeric(c(a,b))
		}
		else if (iniMethod == "lm") {
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
