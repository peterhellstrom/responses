# Calculate number of killed prey (NP):
# n = number of individuals in each predator class (adults, juveniles etc.)
# Example for a raptor pair: c(adults, brood size), e.g. c(2,4)
# der = daily energy requirement in g
# pp = proportion of prey type in the diet (usually expressed in terms of biomass!)
# mmp = mean mass of prey type (in grams)
# T = Time, duration of study in days (defaults to 1)
# This function can use single numbers as inputs, or vectors (all vector must be of the same length)

#' @export
NP <- function(n, der, pp, mmp, T=1) {
	inp <- list(der, pp, n)
	crit <- unique(sapply(inp, length))
	if (length(crit) > 1) stop("Input vectors der, pp & n must all be of equal length")
	T * (sum(n * der * pp) /  mmp)
}

# Daily energy requirement
# cast.rate = daily casting rate of pellets
# corr.factor = single value or vector with correction factor (corrects for digested prey, usually based on feeding experiments in captivity)
# mmp = mean mass in grams of prey tyep/species i
# n = number of prey items of type/species i
# N = total number of pellets analyzed
# output = c("sum","ind"). If "sum" only the sum of all prey types is returned. If "ind", the contribution of each prey type i will be returned
#' @export
der <- function(cast.rate=1.1, corr.factor=1, mmp, n, N, output=c("sum", "ind")) {
	output <- match.arg(output)
	if (output == "sum") {
		sum((cast.rate * corr.factor * mmp * n) / N)
	} else {
		(cast.rate * corr.factor * mmp * n) / N
	}
}

# Examples:
# der(mmp=c(31.50, 48.60, 31.30), n=c(10,15,14), N=20)
# Use different correction factors: input a vector of same length as n
# der(mmp=c(31.50, 48.60, 31.30), corr.factor=c(1.05,1.52,1.7), n=c(10,15,14), N=20)

# Predicted daily food requirements, from scaling theory
# After Nagy 1987 + Bozinovic & Medel 1988
# m = body weight (mass) of animal (in grams)
# eff = assimilation efficiency, a proportion (between 0 and 1)
# cal = caloric content of 1 g prey

# Default values: eff=0.769 (estimate for raptors), cal=6.65 (caloric content for 1 g of small mammal prey)
#' @export
der.pred <- function(m=1000, eff=0.769, cal=6.65) {
	fmr <- 10.9 * m^0.64
	fmr / (eff * cal)
}

# Example:
# der.pred(m=800, eff=0.6)

# Heavier prey: fewer items are necessary to fulfil daily energy requirements.
# For heavier prey, asymptote of functional response will always be lower, so the relevant null hypotheses
# is not if asymptote differ between species, but rather if the asymptote differes more or less than expected!
# Expected difference is the weight ratio between prey types.

# Illustrated for grey-sided voles and lemmings, for the same proportion in diet:
# NP(pp=0.6089, der=140, mmp=31.50, n=4)
# NP(pp=0.6089, der=140, mmp=48.60, n=4)

# Input vectors
# NP(pp=c(0.5,0.5), der=c(110, 140), mmp=30, n=c(2,4))

# Example for female, male, nestling (females larger than males)
# NP(mmp=31.50, der=c(120,110,140), n=c(1,1,4), pp=c(0.6,0.6,0.75))

# Assume two predator classes/stages.
# If pooled proportion in diet in terms of biomass is p.hat=0.65, but proportions differ between classes.
# Assume prop1=0.8, then prop2 = 2*p.hat - prop1 = 2*0.65 - 0.8 = 0.5
# Compare estimates of NP, pooled vs. unpooled:
# NP(mmp=31.50, der=c(110,140), n=c(2,4), pp=c(0.65,0.65)) # Both classes have the same proportion = 16.10
 #NP(mmp=31.50, der=c(110,140), n=c(2,4), pp=c(0.5,0.8)) # Difference between prey classes = 17.71
# NP(mmp=31.50, der=c(110,140), n=c(2,4), pp=c(0.8,0.5)) # Switch proportion order = 14.47619
# Difference depends on der and n
