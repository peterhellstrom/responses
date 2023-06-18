
########################################################################
# Common detection functions
# hn = half-normal
# hr = hazard rate
# u = uniform
# Should add negative exponential for completeness

hn <- function(x,sigma) exp(-(x^2)/(2*sigma^2))
hr <- function(x,sigma,b) 1 - exp(-(x/sigma)^-b)



########################################################################
# PLOT FUNCTIONS

plot.dist.hn <- function(dat,object) {
# dat = raw data file
# object = stored object fitted with dist.hn
# Draw the detection function and scaled histogram:
# Create a barplot (to compare with distance output with 30 bins)
# breaks: seq(from=0,length=31,by=4.767)

max.dist <- max(dat$PerpDistance, na.rm=T)
intrvl <- 4.767 # Distance default
bins <- round(max.dist/intrvl,0)

x <- seq(from=0,to=max.dist+intrvl,by=intrvl)

bp1 <- hist(dat$PerpDistance[-object$rm], breaks=x, plot=FALSE)

# Create lines
yv1 <- bp1$density/object$f0
yv2 <- sapply(1:length(yv1), function(i) c(0,rep(yv1[i],times=2),0))

xv <- rep(bp1$breaks,each=4)[-c(1:2)]
xv <- xv[-c(length(xv)-1,length(xv))]

# Create the full plot:
windows(12,8)
plot(x=c(0,max.dist),y=c(0,2.5), type="n", xlab="Distance", ylab="Detection function",
font.lab=2, las=1, main=object$species, xaxt="n")

axis(side=1,at=seq(0,max.dist,20))

lines(xv,yv2,col="blue") # adds bars

abline(h=0,lty=3)
abline(v=0,lty=3)
# Add detection function
curve(hn(x=x,sigma=object$sigma),from=0,to=max.dist,col=2,lwd=2,add=T)
}



########################################################################
# DISTANCE ESTIMATION
# Half-normal detection function (key without adjustment series)
# This closed-form estimator can only be used for the half-normal key,
# otherwise numerical solutions to M-L methods are necessary...
# Not exactly equal to DISTANCE-results (but very close)

# Only point estimation so far,
# no interval estimation

# Required input: an object with fields (columns) PerpDistance, ClusterSize
# Source: section 3.3.4 Buckland et al. 2001

dist.hn <- function(dat,w=FALSE) {
# dat = input object with the data

rm.ind <- which(dat$PerpDistance<0 | is.na(dat$Object) | is.na(dat$PerpDistance) | is.na(dat$ClusterSize) | dat$ClusterSize<0)

x <- dat$PerpDistance[-rm.ind]
cs <- dat$ClusterSize[-rm.ind]

if (is.numeric(w)) w <- w
if (w==FALSE) w <- max(dat$PerpDistance, na.rm=T)

g <- function(x,sigma) exp(-(x^2)/(2*sigma^2))

sumdist <- sum(x^2)
nobs <- length(x)
sigma <- sqrt(sumdist/nobs)

gx <- g(x=x, sigma=sigma)
# mu <- integrate(g, sigma=sigma, lower=0, upper=w)$value # gx dx
mu <- integrate(g, sigma=sigma, lower=0, upper=Inf)$value # gx dx
# mu <- sqrt((pi*sigma^2)/2)
fx <- gx / mu

#f0 <- 1 / mu
f0 <- sqrt(2/(pi*sigma^2))
esw <- mu

xv <- seq(from=0,to=w,by=0.1)
gxv <- g(x=xv,sigma=sigma)
fxv <- gxv / mu

windows(8,4)
par(mfrow=c(1,2))
plot(xv, gxv, font.lab=2, las=1, type="l", ylab="g(x)", main="Detection function")
plot(xv, fxv, font.lab=2, las=1, type="l", ylab="f(x)", main="Probability density function")
abline(h=f0, lty=2, col=2)
par(mfrow=c(1,1))

# Effort
l <- sum(as.vector(tapply(dat$LineLength,dat$LineLabel,sum)/tapply(dat$LineLength,dat$LineLabel,length)))
samples <- length(unique(dat$LineLabel))

a <- 2*w*l
p <- mu / w

# Density of clusters
ds <- 1000*(nobs*f0)/(2*l)

# es (expected cluster size)
avg.s <- mean(cs)
# Calculate size biased estimate of expected cluster size
# see section 3.5.4 in Distance book
xv <- hn(x=x, sigma=sigma)
z <- log(cs)
fm <- lm(z ~ xv)
es <- exp(sum(coef(fm)) + var(z)/2)

# Density
d <- 1000*(nobs*f0*es)/(2*l)

# Short output (compare with DISTANCE)
out <- list(
rm = rm.ind,
Effort = l,
Samples = samples,
Width = w,
Observations = nobs,
sigma = sigma,
f0 = f0,
ESW = esw,
p = p,
area = a,
DS = ds,
mean.S = avg.s,
E.S = es,
D = d)

print(out)
}



########################################################################
# VARIANCE ESTIMATION

# Works only when there are no duplicates in sample layers, so for my
# purpose it only works within years.
# Since variation in encounter rate most often is the largest component
# of variation of a density estimate, it might be good to have a function
# that calculates this component.  

# Calculates Var(n) and Var(n/L) [encounter rate] at the Sample (i.e. line) level
# Input is a data frame that must have at least these three fields:
# Territory, LineLength, LineLabel & Object
# The input data should be a data frame with observations

var.n <- function(dat, level="stratum", pooled=TRUE) {

# Estimator changed to the same used by Distance 6.0 (R2)
# assumes that each line is a valid replicate (tried to re-write to allow subsamples etc,
# but that would require a complete re-write of the function...
# Count number of detections in each strata
# Only works if stratum is single-level (e.g. Site), doesn't work for e.g. Species*Site
# You have to filter data manually first if you have this setup.

if (level=="substratum" & is.null(dat$SubStratum)) { stop("No SubStratum layer in this file, exit") }

obs <- data.frame("LineLabel"=dat$LineLabel,"Object"=dat$Object,"PerpDistance"=dat$PerpDistance,"ClusterSize"=dat$ClusterSize)

if (level=="sample") {
strat <- unique(data.frame("Territory"=dat$Territory,"LineLabel"=dat$LineLabel,
"Stratum"=dat$Stratum, "SubStratum"=dat$SubStratum))
strat <- strat[order(strat$LineLabel),] }

if (level=="substratum") {
strat <- unique(data.frame("LineLabel"=dat$LineLabel,"Stratum"=dat$Stratum, "SubStratum"=dat$SubStratum))
strat <- strat[order(strat$LineLabel),] }

if (level=="global" | level=="stratum") {
strat <- unique(data.frame("LineLabel"=dat$LineLabel,"Stratum"=dat$Stratum))
strat <- strat[order(strat$LineLabel),] }

samp <- unique(data.frame("Territory"=dat$Territory))
samp <- samp[order(samp),]

subsamp <- unique(data.frame("Territory"=dat$Territory,"LineLabel"=dat$LineLabel,"LineLength"=dat$LineLength,"Observer"=dat$Observer))
subsamp <- subsamp[order(subsamp$LineLabel),]


nobs.subsamp <- as.table(ftable(obs$Object~obs$LineLabel))
if (colnames(nobs.subsamp)[1]=="") nobs.subsamp <- nobs.subsamp[,c(-1)]

if (ncol(nobs.subsamp)>1) nobs.subsamp <- as.matrix(apply(nobs.subsamp,1,sum)) 
if (ncol(nobs.subsamp)==1) nobs.subsamp <- as.matrix(nobs.subsamp)

len <- subsamp$LineLength
encounter.rate <- nobs.subsamp/subsamp$LineLength

# Global level
if (level=="global") {
L <- sum(subsamp$LineLength)
n <- sum(nobs.subsamp)

k <- nrow(subsamp) # assumes that each line is a valid replicate

var.nL <- as.numeric((k/(L^2*(k))) * sum(len^2*((encounter.rate - n/L)^2)))

var.n <- var.nL * L^2
nstrat <- length(var.n)
}

# Stratum | Substratum level
else if (level=="stratum" | level=="substratum" | level=="sample") {

  if (level=="stratum") z <- strat$Stratum
  else if (level=="substratum") z <- strat$SubStratum
  else if (level=="sample") z <- strat$Territory

L <- tapply(len,z,sum)
k <- table(factor(z))

len <- split(subsamp$LineLength,z)
nstrat <- length(len)

n <- tapply(nobs.subsamp,z,sum)
encounter.rate <- split(encounter.rate,z)
# Estimator R2
var.nL <- sapply(1:nstrat, function(i) {
(k[i]/(L[i]^2*(k[i]-1))) * sum((len[[i]])^2*((encounter.rate[[i]] - n[i]/L[i])^2))
})

var.n <- var.nL * L^2
}

else stop("Only levels available are global, stratum, substratum or sample")

# Calculate statistics: 
se.n <- sqrt(var.n)
cv.n <- as.numeric(100*se.n / n)

# Encounter rate
nL <- as.numeric(n/L)
se.nL <- sqrt(var.nL)
cv.nL <- 100*(se.nL / nL) 
# "Clumping factor"
poisson.ratio <- as.numeric(var.n / n)

# Create output (four matrices that are brought together in a list in the final step)

rowname <- names(var.n)

Effort <- matrix(nrow=nstrat,ncol=2)
colnames(Effort) <- c("L","k")
Effort[,1] <- L; Effort[,2] <- k;

Detections <- matrix(nrow=nstrat,ncol=4)
colnames(Detections) <- c("n","Var(n)","SE(n)","CV(n)")
Detections[,1] <- n; Detections[,2] <- var.n; Detections[,3] <- se.n; Detections[,4] <- cv.n; 

EncounterRate <- matrix(nrow=nstrat,ncol=4)
colnames(EncounterRate) <- c("n/L","Var(n/L)","SE(n/L)","CV(n/L)")
EncounterRate[,1] <- nL; EncounterRate[,2] <- var.nL; EncounterRate[,3] <- se.nL; EncounterRate[,4] <- cv.nL; 

PoissonRatio <- matrix(nrow=nstrat,ncol=1)
colnames(PoissonRatio) <- c("Poisson ratio")
PoissonRatio[,1] <- poisson.ratio

if (level=="stratum" | level=="substratum" | level=="sample") {
rownames(Effort) <- names(L)
rownames(Detections) <- names(L)
rownames(EncounterRate) <- names(L)
rownames(PoissonRatio) <- names(L)
}

# Output list
out <- list(
Level = level,
Effort = Effort,
Detections = Detections,
EncounterRate = EncounterRate,
VarianceRatio = PoissonRatio
)

print(out)
}



########################################################################

# Two functions that extracts data from the stats file

split.line <- function(line) {
# Be careful with this function! Parsing might go wrong.
# It still gives some errors,
# I should probably update the function to store the number of the row
# where parsing fails, so that I could go back and check what went wrong.
# Another option is to edit the file in Excel first, safer but time consuming

# Line is an indexing vector, called from the read.stats.file function

# Was changed slightly (substring values differed)
# takes each line and returns the different values for stratum, samp, estimator, ...
  if(nchar(line)<6) stop ("This isn't a stats file!")
  stratum <- as.integer(substr(line,1,6))
  if(is.na(stratum)) {
    #this means that 
    ok=FALSE
    samp <- NULL
    estimator <- NULL
    module <- NULL
    statistic <- NULL
    value <- NULL  
    cv <- NULL
    lcl <- NULL
    ucl <- NULL
    degrees.freedom <- NULL
  } else {
    ok=TRUE
#       Parse the rest of the line when ok is true
    samp <- as.integer(substr(line, 8, 13))
    estimator <- as.integer(substr(line,14,15))
    module <- as.integer(substr(line, 16,17))
    statistic <- as.integer(substr(line,18,20))
    value <- as.numeric(substr(line, 21,36))
    # Index values were not correct from here:
    cv <- as.numeric(substr(line, 37,52)) # was 37,45
    lcl <- as.numeric(substr(line, 53,67)) # was 46,60
    ucl <- as.numeric(substr(line, 68,82)) # was 61,75
    degrees.freedom <- as.numeric(substr(line, 83,95)) # was 76,91 - also changed from as.integer to as.numeric
   }
  return(list(ok=ok, stratum=stratum, samp=samp,
              estimator=estimator, module=module, statistic=statistic,
              value=value, cv=cv, lcl=lcl, ucl=ucl, degrees.freedom=degrees.freedom))
}
################################################################################

read.stats.file <- function(stat.file.name) {
#   Purpose: extracts results statistics from the MCDS stat file

#Input:
#  stat.file.name - name of file to look in

#Returns list:
#   Dhat.ind    -   density of individuals
#   Dhat.ind.se
#   Dhat.ind.df
#   Es          -   mean cluster size
#   Es.se
#   f0          -   pooled f0
#   f0.se 
#   n           -   number of detected animals
#   n.se
#   L           -   effort
#   nL          -   encounter rate
#   nL.se

#read the file in and test that it has something in it
  lines.v <- readLines(stat.file.name)
  n.lines <- length(lines.v)
  if(n.lines==0) {
    stop ("Nothing to read!")
  }

# go through each line, looking for the results we want and storing them in the appropriate place

  Dhat.ind <- NULL; Dhat.ind.cv <- NULL; Dhat.ind.df <- NULL; Es <- NULL; Es.cv <- NULL; f0 <- NULL;
  f0.cv <- NULL; n <- NULL; n.cv <- NULL; L <- NULL; nL <- NULL; nL.cv <- NULL

  for (line in 1:n.lines) {
    parsed.line <- split.line(lines.v[line])
    if(parsed.line$ok) {
       # Dhat.ind
       if(parsed.line$module==4 & parsed.line$statistic==2){# density of individuals line
          Dhat.ind[parsed.line$stratum] <- parsed.line$value
          Dhat.ind.cv[parsed.line$stratum] <- parsed.line$cv
          Dhat.ind.df[parsed.line$stratum] <- parsed.line$degrees.freedom}
       # Es
       if(parsed.line$module==3 & parsed.line$statistic==4){# size bias adjusted cluster size 
          Es[parsed.line$stratum] <- parsed.line$value
          Es.cv[parsed.line$stratum] <- parsed.line$cv}
       # f0
       if(parsed.line$module==2 & parsed.line$statistic==4){# f0
          f0 <- parsed.line$value
          f0.cv <- parsed.line$cv}
       # n
       if(parsed.line$module==1 & parsed.line$statistic==1){# n
          n[parsed.line$stratum] <- parsed.line$value}
       # L
       if(parsed.line$module==1 & parsed.line$statistic==3){# L
          L[parsed.line$stratum] <- parsed.line$value}
       # nL
       if(parsed.line$module==1 & parsed.line$statistic==4){# nL
          nL[parsed.line$stratum] <- parsed.line$value
          nL.cv[parsed.line$stratum] <- parsed.line$cv}         
    }
  }
  return(list(
          Dhat.ind=Dhat.ind, Dhat.ind.se=Dhat.ind*Dhat.ind.cv, Dhat.ind.df=Dhat.ind.df, 
          Es=Es, Es.se=Es*Es.cv, f0=f0, f0.se=f0*f0.cv,
          n=n, L=L, nL=nL, nL.se=nL*nL.cv
          ))
 }



########################################################################
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


