# COMPARE TWO DENSITY ESTIMATES

# This function calculates the two-sample t-test described in the Distance book (Buckland et al. 2001)
# on p. 84 - 86. Emphasis is on the situation where you haved a pooled f0
# estimate, but it's also possible to do the tests where f0's have been estimated separately.

# In order to do the calculations, you must first have a stats file.
# The stats file can be exported when an analysis is running in Distance -
# Check under Misc Tab in the Model Definition Properties.
# The format of the stats file is described in the program manual
# p. 96-98 & 292-294 (valid for Distance 5.0 release 2).

# The functions split.line and read.stats.file were modified from the
# online supplement to Rexstad's 2007 study on MCDS distance sampling (great bustards).

# The component var(n) is not given by Distance, but can easily be calculated as
# var(n) = var(nL) * L^2 (see section 3.6.2, eq. 3.78, p. 79  and p. 108-109 in Buckland et al. (2001))
# Note that a new estimator of var(n/L) was introduced in Distance 6.0 release 1
################################################################################

# cp2est calculates the t-statistic for two samples

cp2est <- function(x,statistic="t",independent=FALSE,alpha=0.05) {

# x is a list or vector with DISTANCE output results
# statistic is "t" or "z", z should only be used for large sample
# independent: FALSE or TRUE: If FALSE, it is assumed that data were pooled
# for estimation of f0, as I understand it the parameters d, nL & Es 
# should be estimated by stratum if the formulas in section 3.6.5 shall be used.
# alpha is significance level

# Stop if you have more than two estimates
if (length(x$Dhat.ind) > 2) stop("Only two density estimates can be compared")

# Necessary inputs if independent == TRUE
d1 <- x$Dhat.ind[1]
d2 <- x$Dhat.ind[2]
var.d1 <- (x$Dhat.ind.se^2)[1]
var.d2 <- (x$Dhat.ind.se^2)[2]
df1 <- x$Dhat.ind.df[1] 
df2 <- x$Dhat.ind.df[2]

# if independent == FALSE you also need
n1 <- x$n[1]
n2 <- x$n[2]

var.n1 <- (((x$nL.se)^2) * x$L^2)[1]
var.n2 <- (((x$nL.se)^2) * x$L^2)[2]

f0 <- x$f0
var.f0 <- x$f0.se^2

if (length(x$Es) != 2) { stop("You have estimated Es pooled, quitting") }
  E1 <- x$Es[1]
  E2 <- x$Es[2]
  var.E1 <- x$Es.se[1]^2
  var.E2 <- x$Es.se[2]^2 

L1 <- x$L[1]
L2 <- x$L[2]

# STEP 1: Calculate variances using the delta method
if (independent == FALSE) {
# internal functions 
M <- function(n,E,L) { (n*E) / (2*L) } # eq 3.107
var.M <- function(M,n,var.n,E,var.E) { M^2 * ( (var.n/n^2) + (var.E/E^2) ) } # eq 3.11

M1 <- M(n=n1, E=E1, L=L1)
M2 <- M(n=n2, E=E2, L=L2)

var.M1 <- var.M(M=M1, n=n1, var.n=var.n1, E=E1, var.E=var.E1)
var.M2 <- var.M(M=M2, n=n2, var.n=var.n2, E=E2, var.E=var.E2)

M.diff <- M1 - M2
var.M.diff <- var.M1 + var.M2 # eq 3.110

d.diff <- 1000*(M1-M2)*f0 # eq 3.108, units assumed to be square kilometer and distances in metres
var.d.diff <- d.diff^2 * ((var.M.diff/M.diff^2) + (var.f0/f0^2)) # eq 3.109
}

else if (independent == TRUE) {
d.diff <- d1 - d2
var.d.diff <- var.d1 + var.d2
}

else stop("Variable independent can only be TRUE or FALSE")

# STEP 2: Calculate test statistics and output
if (statistic == "t") {
t <- abs(d.diff / sqrt(var.d.diff))
conf.t <- numeric(2)
t.df <- ((var.d1 + var.d2)^2) / ((var.d1^2 / df1) + (var.d2^2 / df2))
t.crit <- qt(p=(1-alpha/2), df=t.df)
t.p <- 2*(1 - pt(q=t, df=t.df))
conf.t[2] <- d.diff + t.crit*sqrt(var.d.diff)
conf.t[1] <- d.diff - t.crit*sqrt(var.d.diff)
list(d=x$Dhat.ind, d.diff=d.diff, var.d.diff=var.d.diff, 
t=t, df=t.df, t.crit=t.crit, p=t.p, CI.diff=conf.t) 
}

else if (statistic == "z") {
z <- abs(d.diff / sqrt(var.d.diff))
conf.z <- numeric(2)
z.crit <- qnorm(p=(1-alpha/2),mean=0,sd=1)
z.p <- 2*(1 - pnorm(q=z, mean=0, sd=1))
conf.z[2] <- d.diff + z.crit*sqrt(var.d.diff)
conf.z[1] <- d.diff - z.crit*sqrt(var.d.diff)
list(d=x$Dhat.ind, d.diff=d.diff, var.d.diff=var.d.diff, z=z, z.crit=z.crit, p=z.p, CI.diff=conf.z)
}

else stop("Only statistic t or z allowed")
  
}
