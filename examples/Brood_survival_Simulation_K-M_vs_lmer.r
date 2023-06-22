library(responses)
library(lme4)

# Code for simulating two "populations" and test for difference in survival
# between nestlings stages and interaction stage*population.

# Set up a loop for storing p-value for test of interaction effect stage*population
# I assume that clutch size follows a normal distribution (rounded values)
# Survivors are drawn from clutch resp. hatchling stage by binomial sampling.

# But initial results are indeed promising regarding point estimates:
# Calculations from mean, lmer, and Kaplan-Meier are congruent, and estimated CI's contain simulated input parameters.

n.sim <- 150
ps <- numeric(n.sim)
out.h1 <- matrix(nrow=n.sim,ncol=4)
out.f1 <- matrix(nrow=n.sim,ncol=4)
out.h2 <- matrix(nrow=n.sim,ncol=4)
out.f2 <- matrix(nrow=n.sim,ncol=4)

colnames(out.h1) <- c("Simulated","Means","LMER","KM")
colnames(out.f1) <- c("Simulated","Means","LMER","KM")
colnames(out.h2) <- c("Simulated","Means","LMER","KM")
colnames(out.f2) <- c("Simulated","Means","LMER","KM")

# Set up simulation parameters
eggs <- rpois(50,4.445)
# time intervals = 5
# interval specific survival estimates, assume equally spaced
surv <- c(0.8, 0.9, 0.9, 0.95, 0.95)
# Survival function
cumprod(surv)
# Simulate survival
rbinom(eggs,eggs,surv)

# POPULATION 1
# Input parameters for clutch size
n1 <- 95
x1 <- 4.43
sd1 <- 1.19
# Survival probabilities
p11 <- 0.8
p12 <- 0.9

# POPULATION 2
# Input parameters for clutch size
n2 <- 35
x2 <- 3.74
sd2 <- 1.24
# Survival probabilites
p21 <- 0.8
p22 <- 0.9

# Create factors
Population <- factor(rep(c(1,2), times=c(n1,n2)))
Brood <- factor(c(1:n1,1:n2))
Brood <- Population:Brood

# Start for-loop
for (i in 1:n.sim) {

cat("Wait... iteration number", i, "\n")

# Population 1
c1 <- round(rnorm(n1,x1,sd1),0) # Assume that clutch size is normal.
# Survival probabilities (p ij, i=population, j=stage [stage 1 = egg-hatch, stage 2 = hatch-fledge)
h1 <- rbinom(c1,c1,p11) # Calculate number of hatchlings (binomial sampling)
f1 <- rbinom(h1,h1,p12) # Calculate number of fledglings (binomial sampling)
pop1 <- data.frame(c1,h1,f1) # Check that it looks OK

# Repeat the same procedure for population 2:
c2 <- round(rnorm(n2,x2,sd2),0)
h2 <- rbinom(c2,c2,p21)
f2 <- rbinom(h2,h2,p22)
pop2 <- data.frame(c2,h2,f2)

# Next step, create data frame for analysis
# Bind the data together
sim.dat <- data.frame(ClutchSize=c(c1,c2),Hatchlings=c(h1,h2),Fledglings=c(f1,f2),Population,Brood)

# Re-format the data (binomial input)
z <- 3
temp <- lapply(2:z, function(i) data.frame(
	Alive = sim.dat[,i],
	Dead = sim.dat[,i-1] - sim.dat[,i],
	Brood = sim.dat$Brood,
	Population = factor(i))
)

temp <- do.call("rbind", temp)

# Reformat the data once again, a list that later is called by lmer
inp <- list(
	surv = cbind(Alive = temp$Alive, Dead = temp$Dead),
	brood = factor(temp$Brood),
	stage = factor(temp$Time),
	year = factor(temp$Year)
)

# Call lmer
rep.bn.lmer <- glmer(surv ~ stage * year + (1|brood), family=binomial, data=inp)
glmer(surv ~ stage + (1|year), family=binomial, data=inp)

coefs <- rep.bn.lmer@fixef

rep.bn.lmer2 <- glmer(surv ~ stage + year + (1|brood), family=binomial, data=inp)
anova(rep.bn.lmer,rep.bn.lmer2) # Should really use REML=FALSE for fixed effects comparison!!!

pval <- anova(rep.bn.lmer,rep.bn.lmer2)[2,7]
ps[i] <- pval # send p-value to a vector, de-select this step if not running the for loop!


# Calculate means of nestlings per breeding stage ----
means <- cbind(
tapply(sim.dat$ClutchSize,sim.dat$Year,mean),
tapply(sim.dat$Hatchlings,sim.dat$Year,mean),
tapply(sim.dat$Fledglings,sim.dat$Year,mean)
)

#t(
#tapply(temp$Alive,list(temp$Time,temp$Year),sum)/
#tapply(temp$Alive+temp$Dead,list(temp$Time,temp$Year),sum))

# Output comparing simulated values, calculated from means, estimated with lmer and estimated with KM.
out <-
data.frame(
simulated = c(p11,p12,p21,p22), # print simulated p's

from.means = c(
sapply(2:3, function(i) means[1,i]/means[1,i-1]),
sapply(2:3, function(i) means[2,i]/means[2,i-1])
), # Calculate survival p's from means

est.lmer = c(
logIt(as.numeric(coefs[1])),
logIt(as.numeric(coefs[1] + coefs[2])),
logIt(as.numeric(coefs[1] + coefs[3])),
logIt(as.numeric((coefs[1] + coefs[2] + coefs[3] + coefs[4])))
), # from lmer

est.KM = c(
surv.prob(pop1)$S.t,
surv.prob(pop2)$S.t)
) # Compare with modified Kaplan-Meier

# Send output to "Big" matrices:
out.h1[i,] <- as.numeric(out[1,])
out.f1[i,] <- as.numeric(out[2,])
out.h2[i,] <- as.numeric(out[3,])
out.f2[i,] <- as.numeric(out[4,])

}

# Analyze simulation output: ----

# Point estimation
h1.est <- apply(out.h1,2,function (x) quantile(x,c(0.025,0.5,0.975),na.rm=T))
f1.est <- apply(out.f1,2,function (x) quantile(x,c(0.025,0.5,0.975),na.rm=T))

h2.est <- apply(out.h2,2,function (x) quantile(x,c(0.025,0.5,0.975),na.rm=T))
f2.est <- apply(out.f2,2,function (x) quantile(x,c(0.025,0.5,0.975),na.rm=T))

h1.est
f1.est

h2.est
f2.est

# Proportion of simulations where interaction effect (stage * year) was significant
length(which(ps<0.05)) / length(which(!is.na(ps)))

# Interval estimation...?
