# Bacterial growth example
# From Statistics for Environmental Engineers (2002, 2nd ed), chapter 46
# by Berthouex & Brown

bact = data.frame(
	D = c(0.042, 0.056, 0.083, 0.167, 0.333, 0.5, 0.667, 1),
	S = c(221, 87, 112, 120, 113, 224, 1569, 2745),
	X = c(1589, 2010, 1993, 1917, 1731, 1787, 676, 122)
)

bact <- bact[-8,] # Data point removed in book example

# Test initial estimates p0 and plot functions:
p0 <- c(0.7, 105, 0.65)
S0 <- 3000

with(bact, {
	plot(D, S, xlim=c(0,1), ylim=c(0,2800))
	points(D, X, pch=16)
})

curve(p0[2]*x / (p0[1] - x), add=TRUE, from=0, to=0.7, col=1)
curve(p0[3] * (S0 - (p0[2]*x / (p0[1] - x))), add=TRUE, from=0, to=0.7, col=2)

# Estimation
fr <- function(p0) {
  
	# Observed responses
	y <- matrix(c(bact$S, bact$X),ncol=2)
	x <- bact$D
	# Predicted responses
	f1 <- p0[2]*x / (p0[1] - x)
	f2 <- p0[3] * (S0 - (p0[2]*x / (p0[1] - x)))
  
	f <- matrix(c(f1,f2), ncol=2)
  
	# Box-Draper criterion (see chapter 4 in Bates & Watts 1988).
	# Minimizes the determinant of the crossproduct of the residual matrix y-f.
	# svd = Singular Value Decomposition of a Matrix
	# prod = Product of Vector Elements
	return(prod(svd(y - f, nu=0, nv=0)$d)^2)
}

# Use nlm to minimize the fr function
# nlm = Non-Linear Minimization
out.nlm <- nlm(fr, p=p0, hessian=TRUE, iterlim=100000)
out.nlm
p.nlm <- out.nlm$estimate

# With constrained optimization nlminb
# Appears very sensitive to upper limits, for instance
# upper = c(.75, 120, .65) give similar estimates to nlm,
# but upper = c(.75, 120, .75) give totally different estimates

out.nlminb <- nlminb(start=p0, objective=fr, lower = c(0,0,0), upper = c(.75,120,.65))
out.nlminb
p.nlminb <- out.nlminb$par

cat("Estimates with nlm:", round(p.nlm,3), "\n")
cat("Estimates with nlminb:", round(p.nlminb,3), "\n")

# Plot curves with estimated parameters
with(bact, {
	plot(D, S, xlim=c(0,1), ylim=c(0,2800))
	points(D, X, pch=16)
})

curve(p.nlm[2]*x / (p.nlm[1] - x), add=TRUE, from=0, to=0.7, col=1, lty=2)
curve(p.nlm[3] * (S0 - (p.nlm[2]*x / (p.nlm[1] - x))), add=TRUE, from=0, to=0.7, col=2, lty=2)

curve(p.nlminb[2]*x / (p.nlminb[1] - x), add=TRUE, from=0, to=0.7, col=1, lty=1)
curve(p.nlminb[3] * (S0 - (p.nlminb[2]*x / (p.nlminb[1] - x))), add=TRUE, from=0, to=0.7, col=2, lty=1)

# nlminb fit appears better!
