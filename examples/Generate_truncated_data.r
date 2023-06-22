################################################################################
source("W:/projects/-= R Development =-/distRibutions/distRibutions.r")
################################################################################
# Additional functions

# This function calls genRandom and replicates the simulation n.repl number of times
# Could have been done directly with replicate...
repl.genRandom <- function (n.repl,distr,para,g,...) replicate(genRandom(distr=distr, para=para, g=g, ...), n=n.steps)

################################################################################

# Generate truncated data
# Check computation time

mu <- 10
sigma <- 3.5

curve(dnorm(x=x, mean=mu, sd=sigma), from=0, to=20, ylim=c(0,0.2), n=1001)
# Equal to
#curve(dtrunc(x=x, spec="norm", mean=mu, sd=sigma), from=0, to=20, col=2, n=1001, add=T)
# Left truncation
a <- 5
curve(dtrunc(x=x, spec="norm", mean=mu, sd=sigma, a=a), from=0, to=20, col=2, add=T, n=1001)

# Test the genRandom function, compare with rnorm and rtrunc
# genRandom is basically a wrapper function for either R-basic r*-functions or rtrunc.

n <- 100000 # sample size

# Only left truncation, b=Inf
system.time(y1 <- rnorm(n, mean=mu, sd=sigma))
system.time(y2 <- rtrunc(n, spec="norm", mean=mu, sd=sigma))
system.time(y3 <- rtrunc(n, spec="norm", a=a, mean=mu, sd=sigma))
system.time(y4 <- genRandom(distr="norm", para=list(n=n, mean=mu, sd=sigma), catgs=1))
system.time(y5 <- genRandom(distr="trunc", para=list(n=n, spec="norm", a=a, mean=mu, sd=sigma), catgs=1))

# Plot the data, check that the functions behave as intended:
dev.new()
plot(x=c(0,20), y=c(0,0.2), type="n", xlab="x", ylab="density")
lines(density(y1),col=1,lty=1)
lines(density(y2),col=2,lty=2)
lines(density(y3),col=3,lty=3)
lines(density(y4),col=4,lty=4)
lines(density(y5),col=5,lty=5)

# Conclusion: 
# rnorm & rtrunc are equally fast
# genRandom is (very) slightly slower for non-truncated data.
# but... genRandom is around 5 times slower for truncated data compared with rtrunc!
# Why?

################################################################################
# Simulate grouped data with a grouping variable

# Group parameter
mu4 <- c(3,4,5,6) # means
sigma4 <- c(1,1,1,1) # sds
a4 <- c(0,0,0,0) # left truncation
b4 <- c(8,8,8,8) # right truncation

# Create grouping variable
g4 <- sample(1:4, n, replace=TRUE)

system.time(y6 <- rtrunc(n, spec="norm", a=a4[g4], b=b4[g4], mean=mu4[g4], sd=sigma4[g4]))
system.time(y7 <- rnorm(n, mean=mu4[g4], sd=sigma4[g4])) # no truncation
system.time(y8 <- genRandom(distr="trunc", para=list(n=n, spec="norm", a=a4, b=b4, mean=mu4, sd=sigma4), g=g4, catgs=4))
# Do not forget the g argument in genRandom!

# Conclusion: rnorm with grouped and ungrouped data: no difference
# rtrunc is 6-8 times slower with a grouping variable than rnorm!
# with a grouping variable, rtrunc & genRandom takes equal time.
y6.split <- split(y6, g4)
y7.split <- split(y7, g4)
y8.split <- split(y8, g4)

sapply(y6.split, range)
sapply(y7.split, range)
sapply(y8.split, range)

sapply(y6.split, mean)
sapply(y7.split, mean)
sapply(y8.split, mean)

# Plot density of data generated with rtrunc, rnorm & genRandom:
dev.new()
plot(x=c(0,10), y=c(0,0.4), type="n", xlab="x", ylab="density", main="rtrunc")
for (i in 1:length(y6.split)) lines(density(y6.split[[i]]), col=i, lty=i)

dev.new()
plot(x=c(0,10), y=c(0,0.4), type="n", xlab="x", ylab="density", main="rnorm")
for (i in 1:length(y7.split)) lines(density(y7.split[[i]]), col=i, lty=i)

dev.new()
plot(x=c(0,10), y=c(0,0.4), type="n", xlab="x", ylab="density", main="genRandom")
for (i in 1:length(y8.split)) lines(density(y8.split[[i]]), col=i, lty=i)

# Histogram for data generated with genRandom:
dev.new()
par(mfrow=c(2,2))
for (i in 1:length(y8.split)) {
	hist(y8.split[[i]], breaks=30, xlab="Power", col="steelblue", freq=FALSE, main=paste("Group",i))
	curve(dtrunc(x=x, spec="norm", mean=mu4[i], sd=sigma4[i], a=a4[i], b=b4[i]), n=1001, from=0, to=10, add=T)
}
par(mfrow=c(1,1))

# with repl.genRandom
n <- n
n.repl <- 6
pv.trunc <- list(n=n, spec="norm", a=a4, b=b4, mean=mu4, sd=sigma4)

# pw.para.trunc <- list(n=n.sites, spec="norm", a=c(0,0,0,0), mean=c(1,2,3,3), sd=c(0.5,0.75,1,1)) # left-truncated only
system.time(y9 <- repl.genRandom(n.repl=n.repl, distr="trunc", para=pv.trunc, g=g4))

################################################################################

