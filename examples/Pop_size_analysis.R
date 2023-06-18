################################################################################
# Source the file containing the necessary functions
setwd("C:/WORK/ANALYSER/-= Estimating abundance =-/AccumulationCurves")
source("AccPopSize.R")
################################################################################
# Last changed December 13 2010
################################################################################

# Analysis of arctic fox data from the Helags area in Sweden, published by Meijer et al. (2008) Animal Conservation.

# 66 scats were genotyped, 30 unique genotypes were found.
# Create a data frame:
samples <- 1:66
genotypes <- rep(1:30, c(3,2,3,7,4,2,3,2,1,3,5,4,2,3,4,1,1,1,1,1,1,1,1,2,3,1,1,1,1,1))
rs <- 1000 # number of resamplings
foxes <- data.frame(cbind(samples,genotypes))
# rm(samples,genotypes)

# Create a text file with the data:
# write.table(foxed,"foxes_rarefaction_data.txt",sep="/t", quote=F, col.names=T, row.names=F)

# Test Chao's (2) estimator:
# The input is capture FREQUENCIES
# The function chao.est can only take single-line vectors as input (not a matrix or data.frame in its current implementation)
genToFreq(genotypes) # Capture frequencies for foxes
chao.est(genToFreq(genotypes)) # Chao 2 (chao.est() takes capture frequencies as input)
chao.est(c(14,5,6,3,1,0,1)) # Empirical estimate

# Bootstrap estimate
# Notice that median estimate is equal to the number of genotypes observed, and that
# the lower confidence limit is much lower than the actual observed number of genotypes.
# It appears that a nonparametric bootstrap is not very good for this data set!?!?
# If sampling with replacement from the original sample, the number of drawn individuals
# is always below the observed number of individuals in the sample.
# Is this because the number of individuals "captured" once is nearly half of the sample?
genToFreq(genotypes)[1] / sum(genToFreq(genotypes))
windows(7,6) # Plot capture frequencies;
barplot(genToFreq(genotypes), names.arg=1:length(genToFreq(genotypes)),
	xlab="Number of recaptures / individual", ylab="Frequency", font.lab=2, las=1,
	main="Capture frequencies")
abline(h=0)

inp.resamp.boot <- resampl.data(samples, genotypes, resamplings=rs, replace=TRUE)
# Calculate number of unique genotypes per resampling with replacement:
nUG.rs <- sapply(1:rs, function(i) length(unique(inp.resamp.boot$y[,i])))
windows(7,6)
hist(nUG.rs, breaks=10, col="lightgrey", xlab="Unique genotypes in resampled data set", font.lab=2, las=1)

n.boot <- t(sapply(1:ncol(inp.resamp.boot$y), function(i) chao.est(genToFreq(inp.resamp.boot$z[,i]))))

hist(n.boot[,1], col="lightgrey", breaks=50, xlab="Estimate", ylab="Frequency", font.lab=2, las=1, main="Bootstrap estimate, Chao2")
abline(v=quantile(n.boot[,1], c(0.025, 0.5, 0.975)), col=c(2,1,2), lty=c(2,1,2), lwd=c(1,2,1))
quantile(n.boot[,1], c(0.025, 0.5, 0.975))

chao.est(c(43,16,8,6,0,2,1)) # Examples from Chao's 1987 paper
chao.est(c(65,12,0,0))
################################################################################
# Sample from a certain distribution?
# What's the prob. that a certain individual is sampled n times?

t1 <- FreqTogen(n=genToFreq(genotypes),freq=1:7)
genToFreq(t1)
genToFreq(genotypes)

mean(genToFreq(genotypes))
x <- rpois(n=length(unique(genotypes)), lambda=mean(genToFreq(genotypes)))
# To be continued... constrain both n and sum(x)?

################################################################################
# Different ways of generating starting values for nls:
# The first step is to evaluate the Kohn function (which is done by the get.ini.mm() function).
# Starting values for Eggert and Chessel are then arbitrarily calculated by the acc.fit() function, as
# a/2 and -b2/1000 (for Eggert), where a = asymptote of Kohn equation and b = half-saturation constant of Kohn equation.
# The function get.ini.mm is an internal function that is called by the main function acc.curve,
# as specified by the iniMethod argument.
# Check compuation times for a single run, the manual input is a shortcut that ignores iterative
# fitting of starting values and is of course the fastest option, followed by lm.
# Best (=most accurate) starting parameters are returned by the option nls.
# These options are mainly important for simulation purposes, where computation time cen be drastically changed
# when using different methods.
system.time(get.ini.mm(samples,cumunique(sample(genotypes)), iniMethod=c(40,40))) # Manual input
system.time(acc.curve(foxes$samples, foxes$genotypes, iniMethod="nls")) # Default, nls Michaelis-Menten method
system.time(acc.curve(foxes$samples, foxes$genotypes, iniMethod="rlm")) # Lineweaver-Burk, robust regression (MASS)
system.time(acc.curve(foxes$samples, foxes$genotypes, iniMethod="lm")) # Lineweaver-Burk, OLS

# Use the acc.curve() function for estimation of a SINGLE sample
# Mainly for learning and ilustrative purposes!
fm1 <- acc.curve(foxes$samples, foxes$genotypes, iniMethod="lm", Log=FALSE)
fm1
fm1$ParameterEstimates # Parameter estimates
fm1$msel # Model selection
# Plot
plot.acc.curve(fm1)
plot.acc.curve(fm1, addPoints=FALSE)

# Expected time in seconds for n resamplings:
system.time(acc.curve(foxes$samples, foxes$genotypes))["elapsed"]*resamplings

################################################################################
# Illustration of the main functions:

# The resampling part with the acc.curve.resamp() function:
# Equal to acc.curve, just also set the number of resamples
# NOTE: the Log() options is included only for exploratory purposes.
# But it is noticeable the Michaelis-Menten model (Kohn) for log-transformed data produces results more similar
# to Eggert and Chessel for untransformed data (i.e. overestimation of N is counteracted by log-transforming when using the M-M equation).

# The resampling and estimation procedure is done in two steps.
# The first step is to generate resampled input data with the function resampl.data():
inp.resamp <- resampl.data(samples, genotypes, rs, Log=FALSE)
str(inp.resamp)

# Some descriptive plots to check that the simulation procedure was succesful:

windows(7,6)
matplot(inp.resamp$y, pch=16, cex=0.6, type="l", col=1, 
xlab="Samples", ylab="Unique genotypes", font.lab=2, las=1, bty="l", main="Resampled data")
points(inp.resamp$x, rowMeans(inp.resamp$y), pch=16, col=2)

windows(7,6)
plot(inp.resamp$x, rowMeans(inp.resamp$y), xlab="Samples", ylab="Unique genotypes", font.lab=2, las=1, pch=16, col=2,
bty="l", main="Resampled data, mean values +/- SD")
arrows(inp.resamp$x, rowMeans(inp.resamp$y) - (apply(inp.resamp$y, 1, sd)/2), 
inp.resamp$x, rowMeans(inp.resamp$y) + (apply(inp.resamp$y, 1, sd)/2), length=0, angle=0)

# Estimation, starting values generated by nls:
fm2 <- acc.curve.resamp(inp.resamp, iniMethod="nls", IC=TRUE, distmom=FALSE)
fm2[1:6]
fm2$ParameterSummary

# Manual input of starting values, the vector iniVals are used for all resampled data sets:
iniVals <- get.ini.mm(inp.resamp$x, rowMeans(inp.resamp$y), iniMethod="nls")
fm3 <- acc.curve.resamp(inp.resamp, iniMethod=iniVals, IC=TRUE, distmom=FALSE)
fm3[1:6]
fm3$ParameterSummary

# SAVE OUTPUT
# dput(fm2, "arctic_fox_analysis")
# list2ascii(fm2[1:8],"af_Helags_1000_resamplings.txt")
################################################################################
# Create various plots: histogram, density and comparison of models:
# Argument model=c("Kohn", "Eggert", "Chessel")

plot.curve(fm2, model="Kohn")
plot.curve(fm2, model="Kohn", ylim=c(0,60))
plot.hist(fm2, model="Kohn", breaks=50)
plot.hist(fm2, model="Kohn", breaks=50, drawN=FALSE)
plot.hist(fm2, model="Kohn", plot.type="density")

plot.curve(fm2, model="Eggert")
plot.hist(fm2, model="Eggert", breaks=50)
plot.hist(fm2, model="Eggert", plot.type="density")

plot.curve(fm2, model="Chessel")
plot.hist(fm2, model="Chessel", breaks=50) # Normality of estimates is almost perfect...
plot.hist(fm2, model="Chessel", plot.type="density")

# Plot a model comparison
plot.est.comp(fm2)

# Check the model with the highest estimate (for Kohn)
modEsts <- fm2$Coefficients[,"a.Kohn"]
inds <- which(modEsts==max(modEsts))
modEsts[inds]

fm.test <- acc.curve(samples=inp.resamp$x, genotypes=inp.resamp$y[,inds], shuffle=FALSE)
fm.test
plot.acc.curve(fm.test)
plot.acc.curve(fm.test,plot.type="extended",xlim=c(0,2000))

################################################################################
# Can the large error with Kohn be caused by an extreme correlation in model parameters?
# And is it possible to account for this correlation in R to get unbiased estimates?
plot.par.cor(fm2)

# Also, can the distribution of parameters
library(fitdistrplus)
library(gamlss)

plotdist(fm2$Coefficients[,"a.Kohn"], breaks=50, col="steelblue")
plot(fitdist(fm2$Coefficients[,"a.Kohn"], "lnorm"), breaks=50, col="steelblue")
plot(fitdist(fm2$Coefficients[,"a.Kohn"], "gamma"), breaks=50, col="steelblue")

histDist(fm2$Coefficients[,"a.Kohn"], family="ST4", nbins=100)

plotdist(fm2$Coefficients[,"a.Eggert"], breaks=50, col="steelblue")
plot(fitdist(fm2$Coefficients[,"a.Eggert"], "lnorm"), breaks=50, col="steelblue")
plot(fitdist(fm2$Coefficients[,"a.Eggert"], "gamma"), breaks=50, col="steelblue")

################################################################################
# Non-linear mixed model:
# system.time(acc.curve.mixed(inp.resamp))
fm4 <- acc.curve.mixed(inp.resamp)
fm4[1:5]

plot.hist(fm4, model="Kohn")
plot.hist(fm4, model="Kohn", plot.type="density")

plot.hist(fm4, model="Eggert")
plot.hist(fm4, model="Eggert", plot.type="density")

plot.hist(fm4, model="Chessel")
plot.hist(fm4, model="Chessel", drawN=FALSE) # Do not overlay normal distribution on histogram
plot.hist(fm4, model="Chessel", plot.type="density")

plot.est.comp(fm4)

# Can variance components and random effects actually tell us something here?
# Solution to overestimation? Weighting?
################################################################################
# RANDOM QUESTIONS:
# Chessel: underestimation? lower CI-limit always very close to observed number of individuals (and even below)
# What if estimate is lower than observed number of genotypes?
# At which x-values is asymptote reached?
# Movements of animals, territoriality and slope of curve?
# Study design might differ for different animal species.
# Random sampling or target sampling at dens.
# EVALUATE FUNCTIONS: which slope of first derivative is necessary to get unbiased estimates?
# Relationship between capture frequencies and slope of the accumulation curve
# Create a function for calculation of capture frequencies?
################################################################################
# INITIAL (STARTING) VALUES:
# Get initial values - comparison of different options:
get.ini.mm(inp.resamp$x, inp.resamp$y[,1])
t(sapply(1:10, function(i) get.ini.mm(inp.resamp$x, inp.resamp$y[,i])))

get.ini.mm(inp.resamp$x, rowMeans(inp.resamp$y))
get.ini.mm(inp.resamp$yGr[,"samples"], inp.resamp$yGr[,"genotypes"])
get.ini.mm(x=inp.resamp$yGr, iniMethod="grouped")

################################################################################
# COMPUTATION TIME and CODE OPTIMIZATION
# 1000 thousand iterations takes ~1 minute, if iniMethod="nls"
# Setting IC == FALSE saves some time, but not much:
system.time(acc.curve.resamp(inp.resamp, iniMethod="nls", IC=FALSE, distmom=FALSE))
system.time(acc.curve.resamp(inp.resamp, iniMethod="nls", IC=TRUE, distmom=FALSE))

# Is it possible to speed this up?
# Using iniMethod="lm" takes ~36 seconds, a considerable improvement in computation time
system.time(acc.curve.resamp(inp.resamp, iniMethod="lm", IC=FALSE, distmom=FALSE))
# The fastest options is to manually input start values, takes ~30 seconds:
iniVals <- get.ini.mm(samples,cumunique(sample(genotypes)), iniMethod="nls")
system.time(acc.curve.resamp(inp.resamp, iniMethod=iniVals, IC=FALSE, distmom=FALSE)) # ~30 s
system.time(acc.curve.resamp(inp.resamp, iniMethod=iniVals, IC=TRUE, distmom=FALSE)) # ~36 s
# This also shows that roughly half the computation time is spent on calculating starting values
# when setting iniMethod="nls".

################################################################################
	