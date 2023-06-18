# Function that obtains likelihood-profiles and confidence intervals for nls- or mle2-objects
pci <- function(object, plot=FALSE, output=c("simple","full")) {
	
	if (class(object) == "lm") stop("object class lm is not implemented")
	else if (class(object) == "glm") stop("object class glm is not supported")

	#p <- profile(object)
	p <- object
	if (plot) plot(profile(object))
	output <- match.arg(output)
	
	if (output == "simple") {
		if (class(object)=="mle2") {
			out <- data.frame(
				'coef' = object@coef,
				'se' = summary(object)@coef[,2],
				'ci' = confint(p)
			)
		} else if (class(object)=="nls") {
			out <- data.frame(
				'coef' = coef(object),
				'se' = summary(object)$coef[,2],
				'ci' = confint(p)
				)
		}
	} else if (output == "full") {
		if (class(object)=="mle2") {
			out <- list(
				'coef' = object@coef,
				'vcov' = vcov(object),
				'cov2cor' = cov2cor(vcov(object)),
				'se' = summary(object)@coef[,2],
				'ci' = t(confint(p))
				)
		} else if (class(object)=="nls") {
			out <- list(
				'coef' = coef(object),
				'vcov' = vcov(object),
				'cov2cor' = cov2cor(vcov(object)),
				'se' = summary(object)$coef[,2],
				'ci' = t(confint(p))
			)
		}
	}
	
	out
}

