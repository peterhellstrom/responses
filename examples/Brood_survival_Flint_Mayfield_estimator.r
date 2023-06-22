# Mayfield-estimator
dat1 <- structure(list(
	X0 = c(4, 4, 2, 4, 4, 4, 5, 4, 2, 5),
	X10 = c(3, 2, 2, 4, 4, 2, 6, 4, 2, 4),
	X20 = c(3, 3, 3, 6, 3, 2, 5, 3, 2, 4),
	X30 = c(3, 2, 3, 6, 4, 2, 4, 2, 2, 3)
	), .Names = c("X0", "X10", "X20", "X30"), 
	row.names = c(NA, -10L), class = "data.frame")


dat2 <- structure(list(X1 = c(NA, NA, NA, 4, NA, 4, NA, NA, NA, NA), 
    X2 = c(NA, NA, NA, 4, 4, 4, NA, NA, 2, NA), X3 = c(4, NA, 
    2, NA, NA, NA, NA, NA, 2, NA), X4 = c(NA, NA, 2, NA, NA, 
    NA, NA, NA, 2, NA), X5 = c(NA, 4, 3, NA, NA, NA, 5, NA, 2, 
    5), X6 = c(NA, 3, NA, NA, NA, NA, NA, NA, NA, NA), X7 = c(NA, 
    2, NA, NA, 4, NA, 6, NA, 2, NA), X8 = c(NA, NA, NA, NA, 4, 
    NA, 6, NA, NA, NA), X9 = c(NA, 2, NA, NA, 4, NA, NA, NA, 
    NA, NA), X10 = c(NA, NA, NA, NA, NA, 2, NA, NA, NA, NA), 
    X11 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), X12 = c(NA, 
    NA, NA, NA, NA, NA, NA, NA, NA, NA), X13 = c(NA, NA, NA, 
    NA, NA, NA, NA, NA, NA, NA), X14 = c(NA, NA, NA, NA, 4, NA, 
    6, 4, NA, 4), X15 = c(3, NA, NA, NA, NA, 2, 5, NA, 2, NA), 
    X16 = c(3, NA, NA, NA, NA, 2, NA, NA, 2, NA), X17 = c(NA, 
    2, NA, NA, 4, 2, NA, NA, NA, NA), X18 = c(3, NA, NA, NA, 
    NA, 2, NA, NA, NA, NA), X19 = c(3, NA, NA, 5, NA, NA, 5, 
    NA, NA, NA), X20 = c(3, 3, NA, NA, NA, NA, 5, NA, NA, NA), 
    X21 = c(NA, NA, NA, 6, NA, NA, 5, NA, NA, NA), X22 = c(NA, 
    NA, NA, 6, NA, NA, NA, NA, NA, NA), X23 = c(3, 2, NA, 6, 
    3, NA, 6, 2, 2, NA), X24 = c(3, 2, 3, 6, NA, 2, 4, NA, NA, 
    NA), X25 = c(3, 2, NA, 6, 4, NA, 4, NA, 2, NA), X26 = c(NA, 
    NA, 4, 6, 4, 2, 4, 2, 2, 3), X27 = c(NA, NA, 3, 6, 4, NA, 
    4, 2, 2, NA)), .Names = c("X1", "X2", "X3", "X4", "X5", "X6", 
"X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", 
"X17", "X18", "X19", "X20", "X21", "X22", "X23", "X24", "X25", 
"X26", "X27"), row.names = c(NA, -10L), class = "data.frame")

library(responses)

mayfield(dat1, method="exposure-days")
mayfield(dat1, method="bart-robson")

mayfield(dat2, method="exposure-days")
mayfield(dat2, method="bart-robson")

#######################################################################
# Plot analysis of dat2
str1 <- mayfield(dat2, method="exposure-days")
xv <- seq(0,50,0.1)
yv <- str1$PointEstimate["dsr"]^xv
S30 <- round(as.numeric(str1$PointEstimate["dsr"]^30),4)
S <- S30

# Confidence interval
cat(paste("Point estimate is", S30, "\n"))
c(round(str1$IntervalEstimate["cil"]^30,4), round(str1$IntervalEstimate["ciu"]^30,4))

plot(xv,yv, xlab="Time (days)", ylab="Brood survival", main=paste("Mayfield estimator, dsr = ",round(str1$PointEstimate["dsr"],4)), font.lab=2, las=1, type="l", lwd=2)
abline(v=0,lty=2,col=2)
abline(v=30,lty=2,col=2)
lines(x=c(-2,30), y=rep(str1$PointEstimate["dsr"]^30,2), lty=3, col=1)

mtext(bquote(S[30] == .(S)),line=0.25)
