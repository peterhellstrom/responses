# Extract coefficients from a gamlss-object.

#' @export
coef.Gamlss <- function(object) {
	family <- object$family[1]
	pars.names <- object$parameters

	pars <- get.coef(object, "coef")
	pars.link <- get.coef(object, "link")
	npar <- length(pars.names); pars.tr <- numeric(npar)

	# Links: identity, log, logit. Not implemented yet: probit, cloglog
	invlogit <- function(x) exp(x)/(1 + exp(x))
	logit <- function(x) log(x / (1-x))

	if ( any(na.omit(pars.link) %in% c("identity","log","logit") == FALSE)) stop("Only 'identity', 'log' and 'logit' link functions are implemented yet")

	for (i in 1:npar) {
		pars.tr[i] <- switch(pars.link[i],
			'identity' = pars[i],
			'log' = log(pars[i]),
			'logit' = logit(pars[i]),
			'NA' = pars[i]) }

	names(pars) <- names(pars.tr) <- names(pars.link) <- pars.names
	data.frame(estimate=pars.tr, transformed=pars, link=pars.link)

}
