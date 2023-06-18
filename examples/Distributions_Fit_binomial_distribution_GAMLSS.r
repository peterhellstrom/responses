setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions")
source("DistributionsFunctions.r")

x <- rbinom(n=5000, size=25, prob=0.3)
y <- cbind(x, 25-x)
fm1 <- gamlss(y~1, family="BI")

coef(fm1)
fitted(fm1,"mu")[1]

coef.gamlss(fm1)

