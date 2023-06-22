library(responses)

# Parameters (Michaelis-Menten parameterization)
a <- 50
b <- 15
theta <- 4
sigma <- 2
# sample size
n <- 50
# x-range
xlims <- c(0,50)

# Simulate data
x <- runif(n,xlims[1],xlims[2])
#y.hat <- TypeIIIb(x=x,a=a,b=b,theta=theta)
y.hat <- TypeIIIb(x=x,a=a,b=b,theta=theta)
y <- rnorm(y.hat,y.hat,sigma)
x <- x[y>0]
y <- y[y>0]
xy = data.frame(x=x, y=y)
with(xy, plot(x,y))

# Fit models, get initial valus

start.h.mle <- getStart.h(x=x, y=y, theta=theta, method="mle") # Holling, mle
start.h.nls <- getStart.h(x=x, y=y, theta=theta, method="nls") # Holling, nls
start.mm.mle <- getStart.mm(x=x, y=y, method="mle") # Michaelis-Menten
start.mm.nls <- getStart.mm(x=x, y=y, method="nls") # For nls, otherwise [method="mle"] start values also includes sigma
# start.mm <- getStart.mm2(x=x, y=y) # possible to use the argument theta as well, e.g. getStart.mm2(x=x, y=y, theta=3)
start.mm
plot(start.mm)
plot(start.h)
# Manually adjust start values
# start.h$ini.IIIb[4] <- 2

# Fit models
m1 <- with(xy, fr.fit(x, y, method="nls", eq="M-M", nls.start=start.mm.nls))
m2 <- with(xy, fr.fit(x, y, method="mle", eq="M-M", nls.start=start.mm.mle))

m3 <- with(xy, fr.fit(x, y, method="nls", eq="Holling", nls.start=start.h.nls))
# This needs more fine-tuning, very sensitive to start values (generation of start values slightly off)
m4 <- with(xy, fr.fit(x, y, method="mle", eq="Holling", nls.start=start.h.mle))

plot(m1)
plot(m2)
plot(m3)
plot(m4)
	
m1$coefs
m2$coefs
m3$coefs
m4$coefs

aicc(m2$models$fm.IIIb)
coef(m2$models$fm.IIIb)
logLik(m2$models$fm.IIIb)
AIC(m2$models$fm.IIIb)

aicc.n(m1$models); aicc.n(m2$models); aicc.n(m3$models); aicc.n(m4$models)

# Fit models without using wrapper function fr.fit
fm.0 <- mle2(minuslogl=LL.null, start=list(a=mean(y), sigma=sd(y)), data=as.list(xy))
fm.I <- mle2(minuslogl=LL.l, start=list(a=coef(lm(y~x-1)), sigma=sd(y)/2), data=as.list(xy))
fm.II <- mle2(minuslogl=LL.mm, start=list(a=a, b=b, sigma=6.62), data=as.list(xy), fixed=list(theta=1))
fm.IIIa <- mle2(minuslogl=LL.mm, start=list(a=a, b=b, sigma=4), data=as.list(xy), fixed=list(theta=2))
fm.IIIb <- mle2(minuslogl=LL.mm, start=list(a=a, b=b, theta=theta, sigma=sigma), data=as.list(xy))

args(LL.mm)
LL.mm(x,y,a=a,b=b,theta=theta,sigma=sigma)
