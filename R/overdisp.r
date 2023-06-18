# Test for overdispersion (works with lme4 and glmmadmb)
# http://glmm.wikidot.com/faq

overdisp <- function(model) {
	# number of variance parameters in an n-by-n variance-covariance matrix
	vpars <- function(m) {
		nrow(m)*(nrow(m)+1)/2
	}
	model.df <- sum(sapply(VarCorr(model),vpars)) + length(fixef(model))
	rdf <- nrow(model.frame(model)) - model.df
	if (class(model) != "glmmadmb") {
		rp <- residuals(model, type="pearson")
	} else {
		rp <- model$residuals
	}
	Pearson.chisq <- sum(rp^2)
	prat <- Pearson.chisq / rdf
	pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
	c(chisq=Pearson.chisq, ratio=prat, rdf=rdf, p=pval)
}

# From Zuur et al. 2013, p. 138 [A Beginner's Guide to GLM and GLMM with R]
overdisp.zuur <- function(model, data) {
	if (class(model) != "glmmadmb") {
		E1 <- resid(model, type="pearson")
	} else {
		E1 <- model$residuals
	}
	N <- nrow(data)
	p <- length(fixef(model)) + 1
	od <- sum(E1^2) / (N - p)
	od
}

# overdisp.zuur(model=m1.admb, data=dat)
