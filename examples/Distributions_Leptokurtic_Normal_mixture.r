# Create leptokurtic distributions
# Hypothesis: A mixture of normal distributions with the same mean and different variances is always leptokurtic
##################################################
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions")
source("DistributionsFunctions.r")
##################################################

n <- 10000

mu <- c(0,0,0)
sigma <- c(1,2,4)

x <- rnorm(n=n, mean=mu, sd=sigma)
descdist(x, boot=999)

f.norm <- fitdist(x, "norm")
f.norm$estimate

dev.new()
plot(f.norm, breaks=30, col="steelblue")

# Fit leptokurtic distribution (T-family) with gamlss
# f.TF <- gamlss(x ~ 1, family="TF")
f.TF <- histDist(x, "TF", density=TRUE, nbins=30)
coef.Gamlss(f.TF)
fpars.TF <- coef.Gamlss(f.TF)$transformed

plot(f.TF)
dev.new(width=12, height=6); plot.Gamlss(f.TF) # Note heavy tails!

plot(density(x))
curve(dTF(x, mu=fpars.TF[1], sigma=fpars.TF[2], nu=fpars.TF[3]), add=TRUE, n=1001, col=2)

dev.new()
plot(density(x), col=2, lwd=2, ylim=c(0,0.4), xlab="x", main="Mixture of normal distributions ==> Leptokurtosis")
for (i in 1:length(mu)) curve(dnorm(x, mean=mu[i], sd=sigma[i]), lty=(i+1), add=T)
curve(dnorm(x, mean=f.norm$estimate[1], sd=f.norm$estimate[2]), col=4, lty=2, lwd=2, add=T)
curve(dTF(x, mu = fpars.TF[1], sigma = fpars.TF[2], nu = fpars.TF[3], log = FALSE), col=3, lwd=2, lty=4, add=T)

legend("topleft", c("Density estimate", "Fitted normal", "Fitted TF"), col=c(2,4,3), lty=c(1,2,4), lwd=c(2,2,2), bty="n")
title(sub="Mixture of distributions with same mean and different variances")
