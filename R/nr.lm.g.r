#' @export
nr.lm.g <- function(x, y, g) {

	if (length(levels(g)) > 2) stop("Current version only allows two levels in factor g")

	# Plot data
	if (min(x) < 0) {xlims <- c(min(x),max(x))}
	if (min(x) >= 0) {xlims <- c(0,max(x))}
	ylims <- c(0,max(y))
	gnum <- as.numeric(g)
	gnum[which(gnum==2)] <- 16

	# Fit the different models
	fm1 <- lm(y ~ x * g) # Varying intercepts and slopes
	fm2 <- lm(y ~ g:x) # Varying slopes, common intercept
	fm3 <- lm(y ~ x + g) # Varying intercepts, common slope
	fm4 <- lm(y ~ g:x - 1) # No intercept, varying slopes
	fm5 <- lm(y ~ x) # No effect of grouping variable, only of covariate

	op <- par(list(mfrow=c(2,3), mar=c(5,4,3,1)))

	plot(x, y, pch=gnum, main="fm1: Varying intercepts and slopes", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=coef(fm1)[1], b=coef(fm1)[2], lty=1)
	abline(a=coef(fm1)[1] + coef(fm1)[3], b=coef(fm1)[2] + coef(fm1)[4], lty=2)

	plot(x, y, pch=gnum, main="fm2: Varying slopes, common intercept", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=coef(fm2)[1], b=coef(fm2)[2], lty=1)
	abline(a=coef(fm2)[1], b=coef(fm2)[3], lty=2)

	plot(x, y, pch=gnum, main="fm3: Varying intercepts, common slope", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=coef(fm3)[1], b=coef(fm3)[2], lty=1)
	abline(a=coef(fm3)[1]+coef(fm3)[3], b=coef(fm3)[2], lty=2)

	plot(x, y, pch=gnum, main="fm4: No intercept, varying slopes", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(a=0, b=coef(fm4)[1], lty=1)
	abline(a=0, b=coef(fm4)[2], lty=2)

	plot(x, y, pch=gnum, main="fm5: No effect of covariate", xlim=xlims, ylim=ylims, font.lab=2, las=1)
	abline(fm5)

	par(mfrow=c(1,1))
	par(op)

	out.nr.lm <- list(fm1,fm2,fm3,fm4,fm5)
	names(out.nr.lm) <- c("fm1","fm2","fm3","fm4","fm5")
	modsel.nr.lm <- ICtab(fm1,fm2,fm3,fm4,fm5,delta=TRUE,sort=TRUE,weights=TRUE,type="AICc",nobs=length(fitted(fm1)))
	r2.nr.lm <- as.numeric(t(sapply(1:length(out.nr.lm), function(i) summary(out.nr.lm[[i]])$r.squared)))

	ind <- attr(modsel.nr.lm,"row.names")[1]
	ind2 <- which(names(out.nr.lm)==ind)

	# Output
	list(
		fm.nr.lm = out.nr.lm,
		aicc.nr.lm = modsel.nr.lm,
		r2.nr.lm = r2.nr.lm,
		best = anova(out.nr.lm[[ind2]],fm5,test="F")
	)
}
