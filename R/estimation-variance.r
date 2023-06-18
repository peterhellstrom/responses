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
