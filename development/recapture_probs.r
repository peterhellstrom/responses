##################################################
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions")
source("DistributionsFunctions.r")
##################################################
ninds <- 50
nsamples <- 200

####
# Same probability
probs <- rep(1/ninds, ninds)

# If recapture probability depends on group structure:
ngroups <- 2
g.probs <- c(0.5,0.5)
g <- factor(sample(1:ngroups, size=ninds, prob=g.probs, replace=TRUE))
# For simplicity, assume that group 2 has twice the prob. of being detected
probs <- as.numeric(g) / sum(as.numeric(g))

# If recapture is unique for each individual, uniform distribution
probs <- runif(n=ninds, min=0, max=1)
probs <- probs/sum(probs)

# Check
sum(probs)

####
# Generate data
x <- sample(1:ninds, size=nsamples, prob=probs, replace=T)
y <- table(table(x))

# Inclusion probability
incl.prob <- length(unique(x)) / ninds

main.str = paste("Recapture data:", ninds, "individuals in pop.,", nsamples, "samples")
bp <- barplot(y / sum(y), xlab="Number of recaptures", ylab="Density", font.lab=2, las=1, main=main.str)
abline(h=0)
title(sub=paste("Inclusion prob. =", round(incl.prob,3), ", Mean recaptures / ind =", round(mean(table(x)),3)))

E.pois <- dpois(as.numeric(names(y)), lambda=mean(table(x)))
E.tr.pois <- dtrunc(as.numeric(names(y)), spec="pois", a=0, b=Inf, lambda=mean(table(x)))
points(bp, E.pois, col=2, pch=16)
points(bp, E.tr.pois, col=4, pch=15)
legend("topright",c("Zero-truncated Poission","Poisson"), title="Expected", pch=c(15,16), col=c(4,2), cex=0.8, bty="n")

# Add probabilites from a zero-truncated Poisson distribution
# xr <- rtrunc(n=1000, spec="pois", a=0, b=Inf, lambda=mean(table(x)))
# range(xr)
# E.tr.pois <- dtrunc(as.numeric(names(y)), spec="pois", a=0, b=Inf, lambda=mean(table(x)))

# Quite possible with overdispersion, try adding negative binomial
# Can the amount of expected overdispersion be calculated from ninds & nsamples?
# Or is there a relationship between inclusion probability and overdispersion?