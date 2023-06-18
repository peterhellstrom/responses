################################################################################
# Source the file containing the necessary functions
setwd("C:/WORK/ANALYSER/-= Estimating abundance =-/AccumulationCurves")
source("AccPopSize.R")
################################################################################
# Last changed December 13 2010
################################################################################

# Parameters:
# genotypes contains (true) population size, must be a list (see below)
# samples is a vector with sample sizes (i.e. number of sampled scats)
# resamplings = number of different data sets created for a given combination of true population size and sample size

################################################################################
# So far, I've done four different simulations:

# Simulation 1, 2 & 3
# These parameters are relevant for the arctic fox in Sweden

# 1
genotypes <- c(10, 25, 50, 75, 100, 125, 150)
samples <- seq(10, 200, 10)
resamplings <- 100

# Start of a rewrite for the whole simulation process...

# GENERATE DATA
# Simulation definitions
sims <- expand.grid(samples=samples, genotypes=genotypes)
runs <- nrow(sims)
total.runs <- runs*resamplings
sims <- data.frame(sims, resamplings=rep(resamplings, runs))
total.runs

# Resample data
sim.data <- lapply(1:runs, function(i) {
		tsmp <- sims[i,"samples"]
		trsmp <- sims[i,"resamplings"]
		tgnty <- 1:sims[i,"genotypes"]
		resampl.data(samples=1:tsmp, genotypes=tgnty, resamplings=trsmp, replace=TRUE)
		})

sim.data.names <- sapply(1:runs, function(i) paste("s",sims[i,"samples"],"g",sims[i,"genotypes"], sep=""))
names(sim.data) <- sim.data.names

# FIT MODELS
fit.data <- lapply(1:runs, function(i) acc.curve.resamp(input.data=out[[i]], iniMethod="lm", IC=FALSE))
names(fit.data) <- sim.data.names

# This looks very strange...or not...
# Resampling distributions could potentially differ much from biologically realistic
# recapture frequencies.
# Sampling with replacement yields nearly "normal" capture frequencies (Poisson process?) for large samples),
# whereas I would believe that most observed data sets are highly skewed (possibly negative binomial)
# Thus, there's a need for a function that generates the process, and also a separate observation model.
# Bayes... and/or compare distribution of nonparametric resamples with parametric distributions.
j <- 65
sims[j,]
tmp <- apply(sim.data[[j]]$z,2,function(x) genToFreq(x))
tmp[[j]]
tmp.chao <- t(sapply(1:resamplings, function(i) chao.est(tmp[[i]])))
windows(8,6)
hist(tmp.chao[,"n"], breaks=50, freq=F, font.lab=2, las=1,
	col="lightgray", xlab="Estimate", ylab="Density", 
	main="Population size estimate, Chao")
abline(v=quantile(tmp.chao[,"n"], c(0.025, 0.5, 0.975)), col=c(2,1,2), lty=c(2,1,2), lwd=c(1,2,1))

genToFreq(sim.data[[j]]$z[,j])
barplot(genToFreq(sim.data[[j]]$z[,j]))

skew(tmp.chao[,"n"])
kurtosis(tmp.chao[,"n"])
plot(ecdf(tmp.chao[,"n"]))

qqnorm(tmp.chao[,"n"])
qqline(tmp.chao[,"n"])

# Things to investigate:
# Relationship between ratio True pop size / number of samples and skewness/kurtosis of distribution of estimates?
# Or is this determined only by number of samples?

# Interaction between model and sample size (potentially also true pop size) in relation to bias:
# Like: gam(bias ~ model*samplesize*truepop)

# HOw to extract summary statistics for modelling?
# For instance CI of estimates - two options:
# 1) median of each resample i, calculate quantile intervals from medians of n resamples
#) 2) Calculate quantile intervals for each resample, and then work with mean or median on each set of n intervals?
# If n is large, should/is there any difference?

# The critical question: what is the process model that should be used to generate this data???
# No paper in molecular ecology has investigated this!

################################################################################
# 2
genotypes <- list(1:25,1:50,1:75,1:100)
samples <- seq(20,150,10)
resamplings <- 100

# 3
genotypes <- list(1:10,1:25,1:50,1:75,1:100,1:125,1:150)
samples <- seq(20,200,10)
resamplings <- 100


# Simulation 4
# Another scenario, with much larger population sizes
# It takes several hours to run this one...
genotypes <- list(1:100, 1:250, 1:500, 1:1000)
samples <- seq(50, 2000, 50)
resamplings <- 100

################################################################################

# I have saved output and did also extract both mean and median values
fm <- acc.curve.sim(genotypes,samples,resamplings)
fm.avg <- acc.curve.xtr.avg(genotypes,samples,fm,statistic="mean")
# or
fm.avg <- acc.curve.xtr.avg(genotypes,samples,fm,statistic="median")

################################################################################
# Plots (save as postscript or pdf)

acc.sim.avg.plot(genotypes,samples,fm.avg,"kohn","point")
acc.sim.avg.plot(genotypes,samples,fm.avg,"eggert","point")
acc.sim.avg.plot(genotypes,samples,fm.avg,"chessel","point")

acc.sim.bias.plot(genotypes,samples,fm.avg,"kohn")
acc.sim.bias.plot(genotypes,samples,fm.avg,"eggert")
acc.sim.bias.plot(genotypes,samples,fm.avg,"chessel")

acc.sim.bias.ratio.plot(genotypes,samples,fm.avg,"kohn")
acc.sim.bias.ratio.plot(genotypes,samples,fm.avg,"eggert")
acc.sim.bias.ratio.plot(genotypes,samples,fm.avg,"chessel")

################################################################################
################################################################################
# Save simulated data
list2ascii(fm,"ss_10-200_tp_10-150.txt")
list2ascii(fm.avg,"ss_10-200_tp_10-150.txt")

# Change to desired output name
dput(fm,"ss_20-150_tp_25-100")
dput(fm.avg,"ss_20-150_tp_25-100_avg") # Store R-object externally

# Change input name (see below for available objects)
fm <- dget("ss_20-150_tp_25-100") # Read stored R-object

fm.avg <- dget("ss_20-150_tp_25-100_avg")

# Naming of output files
# ss_xx-yy_tp_xz-yz
# The first part, ss, is range of sample sizes
# Second part, tp, is true population size
# I have investigated the following

ss_10-200_tp_10-150
ss_20-150_tp_25-100
ss_20-200_tp_10-150 
ss_50-2000_tp_100-1000

# BE CAREFUL when you read a stored object and store it in R as fm
# (see above). When you extract the average values from the fm object, you
# have to set the variables genotypes, samples, and resamplings to the same
# values the fm object was created with!!! Otherwise --> wrong calculations
# and error messages...   

################################################################################
