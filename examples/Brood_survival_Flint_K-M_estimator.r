require(boot)
source("W:/projects/R/nResponse/nResponse.r")

# Data used by Flint et al:
# Example 1
x <- matrix(ncol=4,nrow=10,byrow=TRUE,
c(4,3,3,3,
4,2,3,2,
2,2,3,3,
4,4,6,6,
4,4,3,4,
4,2,2,2,
5,6,5,4,
4,4,3,2,
2,2,2,2,
5,4,4,3),
dimnames = list('Brood' = seq(1,10,1), 'Time' = seq(0,30,10)))

# Example 2
# matrix without brood switching
x <- matrix(ncol=4,nrow=10,byrow=TRUE,
c(4,3,3,3,
4,2,2,2,
2,2,2,2,
4,4,3,3,
4,4,4,3,
4,2,2,2,
5,5,5,4,
4,4,3,2,
2,2,2,2,
5,4,4,3),
dimnames = list('Brood' = seq(1,10,1), 'Time' = seq(0,30,10)))

# Example 3 
# matrix with staggered entry
x <- matrix(ncol=8, nrow=6, byrow=TRUE,
c(5,4,5,5,4,3,3,3,
NA,NA,NA,NA,3,2,2,1,
4,4,4,4,NA,4,3,3,
NA,NA,NA,3,3,2,1,1,
NA,5,4,4,4,4,4,4,
NA,NA,NA,NA,NA,5,5,5),
dimnames = list('Brood' = seq(1,6,1), 'Time' = seq(0,70,10)))


# Test the "new" functions (2013):
flint.km.nc(x)
flint.km.nc.boot(x)
flint.km(x, nboot=2000)

# Test the "old" functions (2009):
m1 <- surv.prob(x) # point estimates
m1
plot(m1)

test <- surv.prob.boot(x=x, nboot=999) # bootstrap
str(test)
test$SE.S.t
test$SE.S.f

# Bootstrap with a vector...
z <- factor(rep(c(1,2), each=5))
testz <- surv.prob.boot.n(x=x,n.boot=999,vec=z)

sapply(1:length(testz), function(i) t(quantile(testz[[i]]$S.t[,1], c(0.025,0.5,0.975))))
sapply(1:length(testz), function(i) t(quantile(testz[[i]]$S.t[,2], c(0.025,0.5,0.975))))
sapply(1:length(testz), function(i) t(quantile(testz[[i]]$S.t[,3], c(0.025,0.5,0.975))))
