source("C:/WORK/ANALYSER/-= R Development =-/distRibutions/distRibutions.r")
devtools::load_all("W:/PROJEKT/R/responses")
library(lme4)
library(boot)

################################################################################
# Randomization test
n.broods <- 100
N.bar <- 4
#N0 <- rtrunc(n.broods, "pois", a = 2, b = 6, lambda=N.bar)
N0 <- rep(5, n.broods)
surv1 <- c(0.8, 0.9, 0.9, 0.95, 0.95)
surv2 <- c(0.7, 0.7, 0.8, 0.8, 0.9)
surv2 <- surv1
x1 <- ProjectSurvival(N0=N0, T=length(surv1), st=surv1)
x2 <- ProjectSurvival(N0=N0, T=length(surv2), st=surv2)
x12 <- rbind(x1, x2)
g <- as.factor(rep(c("p1","p2"), c(n.broods, n.broods)))

surv.prob(x1)$S.t; surv.prob(x2)$S.t
t(sapply(surv.prob.n(x12,g), "[[", "S.t"))
t(sapply(surv.prob.n(x12,g), "[[", "SE.t"))
d2.obs(x12,g)

surv.prob.2rand(x=x12, vec=g, R=500)

surv.prob.boot(x1)

flint.km(x1)
################################################################################
# ToDO
# Include heterogeneity in survival between broods
# Simulate random effects for broods and intervals

n.broods <- 5000
N.bar <- 4
N0 <- rtrunc(n.broods, "pois", a = 3, b = 6, lambda=N.bar)

surv <- c(0.8, 0.9, 0.9, 0.95, 0.95)
pred <- c(mean(N0), mean(N0)*cumprod(surv))
x <- ProjectSurvival(N0=N0, T=length(surv), st=surv)
rbind(colMeans(x),pred)

# Compare GLM with Flint's modified Kaplan-Meier estimator
nsim <- 100
out <- matrix(ncol=2, nrow=nsim)

test <- sapply(1:nsim, function(i) {

N0 <- rtrunc(n.broods, "pois", a = 3, b = 6, lambda=N.bar)
x1 <- ProjectSurvival(N0=N0, T=length(surv), st=surv)
colnames(x1) <- paste("t",1:6,sep="")

plot(c(1,ncol(x1)), c(0,max(N0)), type="n", main=i)
for (j in 1:nrow(x1)) points(1:ncol(x1), x1[j,], type="l", lty=2)
points(1:ncol(x1), pred, type="p", pch=16, col=2)

y <- cbind('alive' = x1[,"t6"], 'dead' = x1[,"t1"] - x1[,"t6"])
m1 <- glm(y ~ 1, family=binomial(link="logit"), weights=log(x1[,1]))
m2 <- flint.km.nc(x1[,c(1,6)])

coefs <- as.vector(c(expit(coef(m1)), m2["S.t",]))
coefs
})
test <- t(test)

hist(apply(test,1,diff), col="steelblue")


################################################################
# MIXED MODEL with glmer
x <- ProjectSurvival(N0=N0, T=length(surv), st=surv)
# x <- x[-which(apply(x,1,prod)==0),]

z <- lapply(1:(ncol(x)-1), function(i) cbind(x[,i],x[,i+1]))
names(z) <- paste("t",1:length(z),sep="")
y <- do.call(rbind, z)
y.w <- as.vector(y[,1])
y <- cbind('alive'=y[,2], 'dead'=y[,1]-y[,2])
brood <- as.factor(rownames(y))
rownames(y) <- NULL
interval <- as.factor(rep(names(z), each=nrow(x)))
dat <- list(y=y, brood=brood, interval=interval, w=y.w)

m2 <- glmer(y ~ interval-1 + (1|brood) + (1|interval), family=binomial(link="logit"), data=dat)
summary(m2)
coef(m2)
expit(fixef(m2))
expit(confint(m2))
flint.km.nc(x)

m3 <- glm(y ~ interval-1, family=binomial(link="logit"), data=dat)
summary(m3)
expit(coef(m3))
expit(confint(m3))

library(car)

m3.boot <- Boot(m3, R=1999)
summary(m3.boot, high.moments=TRUE)
confint(m3.boot, level=.95, type="perc")
expit(confint(m3.boot, level=.95, type="perc"))
hist(m3.boot, legend="separate")
################################################################
