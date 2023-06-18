##################################################
setwd("C:/WORK/ANALYSER/-= Statistical methods in R =-/Distributions")
source("DistributionsFunctions.r")
##################################################

# Plot truncated normal distribution(s)

a <- c(-1, -1, -1, 0, -2) # lower boundary
b <- c(1, 1, 2, 3, 1) # upper boundary
mu <- c(0, 0.5, 1.5, 2, 0) # mean
sigma <- c(1, 0.6, 0.8, 0.8, 1.5) # standard deviation
cols <- c("red","blue","green","yellow","purple") # colors to use in plot

dev.new(width=7, height=6)
# Set up plot window
plot(x=c(-3,4), y=c(0,1), xlab="x", ylab="Density", font.lab=2, las=1, main="Truncated normal distribution", type="n")
# Draw curves
for (i in 1:length(a)) curve(dtrunc(x, spec="norm", a=a[i], b=b[i], mean=mu[i], sd=sigma[i]), from=a[i], to=b[i], add=T, n=10001, col=cols[i])
# Draw vertical lines for distribution truncation points
for (i in 1:length(a)) {
	xv <- c(a[i],b[i]); yv <- dtrunc(x=xv, spec="norm", a=a[i], b=b[i], mean=mu[i], sd=sigma[i])
	for (j in 1:length(xv)) lines(rep(xv[j],2),c(0,yv[j]),col=cols[i], lty=2)
}
# Create legend with bquote & paste
expr <- lapply(1:length(a), function(i) 
			bquote(paste(a==.(a[i]), ", ", b==.(b[i]), ", ", mu==.(mu[i]), ", ", sigma==.(sigma[i]))))
legend("topleft", do.call("expression", expr), lty=1, col=cols, cex=0.9, bty = "n") 

