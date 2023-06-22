################################################################################
# Model isocline numerical response with the nr-function
library(responses)

# Note that function doesn't work if:
# x and y values contain zero OR if x and y values contain 1 (log ratios are undefined in step 2)

# Control number of data points/years
nsim <- 25

# Simulate rodent dynamics:
# Test different dynamics, first-order, second-order with varying cycle length, etc.
# AR(1) process
x <- arima.sim(n=nsim, list(ar=c(-0.1)), sd=sqrt(0.1)) + 2.1
x <- as.numeric(exp(x)) # back-transform from log-scale
# AR(2) process, 4 year cycle with mean 2.1 and var=0.3
x1 <- arima.sim(n=nsim, list(ar=c(0,-0.7)), sd=sqrt(0.3)) + 2.1
x1 <- as.numeric(exp(x1)) # back-transform from log-scale

x2 <- arima.sim(n=nsim, list(ar=c(0,-0.7)), sd=sqrt(0.5)) + 3.3
x2 <- as.numeric(exp(x2)) # back-transform from log-scale

plot(ts(cbind(x1,x2)))

# Isocline numerical response, type II with normal error added, SD
a <- 5; b <- 10
y1 <- abs(rnorm(length(x),(a*x1)/(b+x1),0.4))
y2 <- abs(rnorm(length(x),(0.6*a*x2)/(3*b+x2),0.4))

#gr(y1,method="r",plot=TRUE)
#gr(y2,method="r",plot=TRUE)

dev.new(width=10, height=5)
par(mfrow=c(1,2))
plot(x1, y1, xlab="Prey", ylab="Predator", xlim=c(0, max(c(x1,x2))), ylim=c(0,max(c(y1,y2))), pch=1)
points(x2,y2,pch=16)
plot(scale(x1), y1, xlab="Prey", ylab="Predator", xlim=range(c(scale(x1),scale(x2))), ylim=c(0,max(c(y1,y2))), pch=1)
points(scale(x2),y2,pch=16)
#curve(5.16*(x+1.37)/(1.37 + (x+1.37)), add=TRUE, col=2)
par(mfrow=c(1,1))

tmp1 <- data.frame(x=c(x1,x2), y=c(y1,y2), g=factor(rep(1:2,each=nsim)))
tmp2 <- data.frame(x=c(scale(x1),scale(x2)), y=c(y1,y2), g=factor(rep(1:2,each=nsim)))

f1 <- gnls(y ~ SSmicmen(x,a,b), start=list(a=c(3,3), b=c(10,10)), params = list(a~g-1, b~g-1), weights = varIdent(form = ~1|g), data=tmp1)
f2 <- gnls(y ~ a*(x-d) / ((b-d) + (x-d)), start=list(a=c(3,3), b=c(0.4,0.4), d=c(-1,-1)), params = list(a~g-1, b~g-1, d~g-1), weights = varIdent(form = ~1|g), data=tmp2)
coef(f1)
coef(f2)

test1 <- nr(x=x1,y=y1,plot=TRUE,method="direct")
test[2]


test <- nr(x=x,y=y,plot=TRUE,method="both")
test[3:4]

# Check output
test$fm.dnr$fm6
summary(test$fm.dnr$fm6)
lm.plot(test$fm.dnr$fm2)
lm.plot(test$fm.dnr$fm3)

summary(test$fm.inr$fm11)
lm.plot(test$fm.inr$fm11)

################################################################################
# With a covariate (two levels)

# Control number of data points/years
nsim <- 15

x1 <- arima.sim(n=nsim, list(ar=c(0,-0.7)), sd=sqrt(0.3)) + 2.1
x1 <- as.numeric(exp(x1)) # back-transform from log-scale

x2 <- arima.sim(n=nsim, list(ar=c(0,-0.7)), sd=sqrt(0.3)) + 2.1
x2 <- as.numeric(exp(x2)) # back-transform from log-scale

# Isocline numerical response, type II with normal error added, SD
a <- 5; b <- 10
y1 <- abs(rnorm(length(x1),(a*x1)/(b+x1),0.4))
y1 <- abs(rnorm(length(x1),(0.5 + 0.1*x1),0.4))

a <- 3; b <- 10
y2 <- abs(rnorm(length(x2),(0.075*x2),0.4))

decade <- factor(rep(c("1970s","2000s"),each=nsim))
g <- decade
y <- c(y1,y2)
x <- c(x1,x2)

nr.lm.g(x=x,y=y,g=g)
################################################################################
