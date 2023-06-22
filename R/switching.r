# Make a function that illustrates switching curves
# c = preference
# b = how preference changes with N1/N2, (b = 2, linear change)

#' @export
switching <- function(c, b, type=c("proportions","ratios"), round=3, x.max=NULL, ...) {

  b.text = round(b,round)
  c.text = round(c,round)

  type <- match.arg(type)
  if (type == "proportions") {
    # Basic plot
    plot(x=c(0,1),y=c(0,1),type="n",
         xlab="Proportion available", ylab="Proportion in diet",
         font.lab=2,bty="l",las=1, ...)
    # Null hypothesis
    curve(switch.curve(x,b=1,c), n=1001, lty=5, lwd=2, col=2, add=TRUE)
    # Switching curve
    curve(switch.curve(x,b,c), n=1001, lty=1, lwd=2, col=1, add=TRUE)
    mtext(paste("Switching, c =",c.text,", b =",b.text))
    legend("bottomright",legend=c("Switching","Null case"),lty=c(1,5),lwd=c(2,2),col=c(1,2),bty="n",cex=1)
  }
  else if (type == "ratios") {
    # Basic plot
    y <- switch.curve.ratios(seq(0,x.max,0.1),b,c)
    plot(x=c(0,x.max),y=c(0,max(y)),type="n",
         xlab="Ratio available",ylab="Ratio in diet",
         font.lab=2,bty="l",las=1,ylim=c(0,max(y)), ...)
    # Null hypothesis
    curve(switch.curve.ratios(x,b=1,c), n=1001, lty=5, lwd=2, col=2, add=TRUE)
    # Switching curve
    curve(switch.curve.ratios(x,b,c), n=1001, lty=1, lwd=2, col=1, add=TRUE)
    mtext(paste("Switching, c =",c.text,", b =",b.text))
    legend("bottomright",legend=c("Switching","Null case"),lty=c(1,5),lwd=c(2,2),col=c(1,2),bty="n",cex=1)
  }
}

#' @export
switch.curve <- function(x, b, c) {
  c*x^b / ((1 - x)^b + c*x^b)
}

# x = proportion available
# c = preference
# b = how preference changes with N1/N2, (b = 2, linear change)

#' @export
switch.curve.ratios <- function(x, b, c) {
  c*x^b
}

