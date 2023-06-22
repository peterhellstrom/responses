# PLOT FUNCTIONS

#' @export
plot.dist.hn <- function(dat, object) {
  # dat = raw data file
  # object = stored object fitted with dist.hn
  # Draw the detection function and scaled histogram:
  # Create a barplot (to compare with distance output with 30 bins)
  # breaks: seq(from=0,length=31,by=4.767)

  max.dist <- max(dat$PerpDistance, na.rm=T)
  intrvl <- 4.767 # Distance default
  bins <- round(max.dist/intrvl,0)

  x <- seq(from=0,to=max.dist+intrvl,by=intrvl)

  bp1 <- hist(dat$PerpDistance[-object$rm], breaks=x, plot=FALSE)

  # Create lines
  yv1 <- bp1$density/object$f0
  yv2 <- sapply(1:length(yv1), function(i) c(0,rep(yv1[i],times=2),0))

  xv <- rep(bp1$breaks,each=4)[-c(1:2)]
  xv <- xv[-c(length(xv)-1,length(xv))]

  # Create the full plot:
  plot(x=c(0,max.dist),y=c(0,2.5), type="n", xlab="Distance", ylab="Detection function",
       font.lab=2, las=1, main=object$species, xaxt="n")

  axis(side=1,at=seq(0,max.dist,20))

  lines(xv,yv2,col="blue") # adds bars

  abline(h=0,lty=3)
  abline(v=0,lty=3)
  # Add detection function
  curve(hn(x=x,sigma=object$sigma),from=0,to=max.dist,col=2,lwd=2,add=T)
}
