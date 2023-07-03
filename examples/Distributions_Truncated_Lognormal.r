# Right-truncation with dtrunc ----

# DO NOT FORGET to enter the a argument, with value a=-Inf.
# Otherwise things get wrong...


n <- 10000
mu <- 0.5
sigma <- 1

x <- rlnorm(n, meanlog = mu, sdlog = sigma)
plot(density(x))
plot(density(log(x)))

curve(dlnorm(x, mu, sigma), n=1001, xlim=c(0,20))
# Right truncation with dtrunc: must set a to zero for lognormal
curve(dtrunc(x, spec="lnorm", a=-Inf, b=10, mu, sigma),add=T,col=2,n=1001) # RIGHT
curve(dtrunc(x, spec="lnorm", b=10, mu, sigma),add=T,col=5,n=1001) # WRONG
curve(dLOGNO(x, mu=mu, sigma=sigma), col=4,add=T, n=1001)
# Right truncated with GAMLSS
gen.trun(10, family="LOGNO", type="right")
curve(dLOGNOtr(x, mu=mu, sigma=sigma), col=3, lty=2, add=T, n=1001, xlim=c(0,10))

# What about normal? ----
n <- 10000
mu <- 7
sigma <- 2

x <- rnorm(n,mean=mu,sd=sigma)
plot(density(x))

curve(dnorm(x, mu, sigma), n=1001, xlim=c(0,20))
# Right truncation with dtrunc: must set a to zero for lognormal
curve(dtrunc(x, spec="norm", a=-Inf, b=10, mu, sigma),add=T,col=2,n=1001)
curve(dNO(x, mu=mu, sigma=sigma), col=4,add=T, n=1001)
# Right truncated with GAMLSS
gen.trun(10, family="NO", type="right")
curve(dNOtr(x, mu=mu, sigma=sigma), col=3, lty=2, add=T, n=1001, to=10)
