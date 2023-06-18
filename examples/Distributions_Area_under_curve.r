source("W:/PROJEKT/R/distributions/distRibutions.r")
source("distRibutions.r")

# Draw density function of a theoretical distribution and add an area under the curve:
auc("norm", para = list(mean = 0, sd = 1), from = -3, to = -1)
auc("norm", para = list(mean = 20, sd = 3), xlim = c(9, 31), from = 20, to = 30)
auc("lnorm", para = list(meanlog = 0.2, sdlog = 0.4), xlim = c(0, 3), from = 1, to = 2) 
auc("exp", para = list(rate = 1/3), xlim = c(0, 8), from =  1, to = 2)

# Find tails of a distribution
qnorm(c(0.025, 0.975), mean = 0, sd = 1)
qnorm(c(0.005, 0.995), mean = 0, sd = 1)
# qnorm(c(0, 1), mean = 0, sd = 1)

auc("norm", para = list(mean = 0, sd = 1), from = qnorm(0.975, mean = 0, sd = 1), to = 6, by = 1/1000)

# Not yet possible to add two areas (i.e. tail probabilites) with the function auc.
# However, auc.poly return a matrix with x and y coordinates for an area. Use this function like this:

curve(dnorm(x, mean = 0, sd = 1), n = 1001, from = -3, to = 3)
polygon(auc.poly(distr = "norm",
		para = list(mean = 0, sd = 1),
		from = -6,
		to = qnorm(0.025, mean = 0, sd = 1),
		by = 1/100), col = "lightgrey", lty = 2) # lower tail
polygon(auc.poly(distr = "norm",
		para = list(mean = 0, sd = 1),
		from = qnorm(0.975, mean = 0, sd = 1),
		to = 6, by = 1/100), col = "lightgrey", lty = 2) # upper tail
lines(x = rep(0,2), y = c(0,dnorm(0,0,1)), col = 2, lty = 2)

