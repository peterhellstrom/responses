devtools::load_all("W:/PROJEKT/R/responses")

# Simulate data:
p0 <- list(a=500, b=10, theta=2, sigma=50, n=25, x.max=50)
with(p0, {
	x <- runif(n, min=0, max=x.max)
	y <- abs((a * x^theta / (b^theta + x^theta))) + rnorm(n, 0, sigma)
	xy <<- data.frame(x,y)
})

dev.new()
with(xy, plot(x,y,bty="l",pch=16))

preview (y~TypeIIIb(x,a,b,theta), xy, start=c(a=500,b=10,theta=2.5), variable = 1)
preview (y~TypeIIIb(x,a,b,theta), xy, start=c(a=500,b=10,theta=2), variable = 1)


# Plot the resulting fit

m1 <- with(xy, fr.start.grid(x,y,"TypeIIIb"))
m1$start # Best estimate from grid search
p1 <- m1$nlsSum$parameters[,1] # re-fitted estimates, the function uses the grid search to re-fit the model
# Parameters from Self-starting function
p2 <- getInitial(y ~ SStypeIIIb(x, a, b, theta), data=xy)

# Compare
p1
p2

m21 <- nls(y ~ TypeIIIb(x,a,b,theta), start=as.list(p1), data=xy)
coef(m21)
m22 <- nls(y ~ TypeIIIb(x,a,b,theta), start=as.list(p2)[1:3], data=xy)
coef(m22)
p2

# nls(y ~ a * x^theta / (b^theta + x^theta), start=list(a=500,b=10,theta=1.6), data=xy)

with(xy, plot(x, y, bty="l", pch=16, xlim=c(0,max(x)), ylim=c(0,max(y)), xlab="Prey density", ylab="Prey killed", main="Simulated Type IIIb-response", font.lab=2))
# Add the start values that minimized sums-of-squares in the grid search
with(as.list(p1), curve(TypeIIIb(x,a,b,theta), add=TRUE, n=1001, type="l", lwd=2, col="red"))
# Add parameter estimates
with(as.list(p2), curve(TypeIIIb(x,a,b,theta), add=TRUE, n=1001, type="l", lwd=2, col="blue"))
# Legend
legend("bottomright", c("Grid search","Self-start"), col=c("red","blue"), lwd=c(2,2), bty="n", title="Generate start-values")
