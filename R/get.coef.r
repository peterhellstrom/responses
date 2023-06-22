# Extract coefficients from a gamlss object
# Check help for link functions for the various parameters (necessary to transform to get the correct estimates)

#' @export
get.coef <- function(object, method="coef") {
	pars.fitted <- object$parameters
	npar <- length(pars.fitted)
	if (method %in% c("coef","link") == FALSE) stop("Method must be 'coef' or 'link'")
	if (method == "coef") str <- paste("as.numeric(fitted(object,'", pars.fitted, "')[1])", sep="")
	if (method == "link") str <- rep(paste("object$", pars.fitted, ".link", sep=""))

	pars <- sapply(1:npar, function(i) {
		val <- eval(parse(text=str[i]))
		ifelse(!is.null(val), val, NA)
		})
	names(pars) <- pars.fitted
	pars
}
