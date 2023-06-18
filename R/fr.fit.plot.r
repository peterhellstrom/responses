# Add a chosen type of functional response to an existing plot:
# (works only for objects fitted with nls, predict is not available for mle2.)
# Object = stored object, type = functional response type

plot.fResponse <- function(object, type=c("Type0","TypeI","TypeII","TypeIIIa","TypeIIIb","all"),
	main=NULL, xlab=NULL, ylab=NULL, lcols=c(1,1,2,3,4), ltys=c(3,1,2,3,4), lwd=2, dev.new=FALSE, ...) {
	
	if (length(lcols) < 5) lcols <- rep(lcols[1],5)
	if (length(ltys) < 5) ltys <- rep(ltys[1],5)
	coefs <- object$coefs
	
	if (missing(type)) {
		type <- "all"
	} else {
		type <- match.arg(type)
	}
	
	if (missing(main)) {
		strMain <- "Functional response"
	} else {
		strMain <- main
	}
	
	if (missing(xlab)) {
		strxlab <- "Prey density"
	} else {
		strxlab <- xlab
	}
	
	if (missing(ylab)) {
		strylab <- "Consumption rate"
	} else {
		strylab <- ylab
	}
	
	plot(object$x, object$y,
		bty="l", font.lab=2,
		xlab=strxlab, ylab=strylab, main=strMain,...)

	if (object$eq == "M-M") {
		if (type=="Type0") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			mtext(type)
		} else if (type=="TypeI") {
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			mtext(type)
		} else if (type=="TypeII") {
			with(as.list(coefs[3,]), curve(TypeII(x, a=a, b=b), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIa") {
			with(as.list(coefs[4,]), curve(TypeIIIa(x, a=a, b=b), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIb") {
			with(as.list(coefs[5,]), curve(TypeIIIb(x, a=a, b=b, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			mtext(type)
		} else if (type=="all") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			with(as.list(coefs[3,]), curve(TypeII(x, a=a, b=b), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			with(as.list(coefs[4,]), curve(TypeIIIa(x, a=a, b=b), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			with(as.list(coefs[5,]), curve(TypeIIIb(x, a=a, b=b, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			legend("bottomright",c("Type0","Type I","Type II","Type IIIa", "Type IIIb"), col=lcols, lty=ltys, lwd=lwd, bty="n")
			mtext("All fitted types")
		}
	} else if (object$eq == "Holling") {
		if (type=="Type0") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			mtext(type)
		} else if (type=="TypeI") {
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			mtext(type)
		} else if (type=="TypeII") {
			with(as.list(coefs[3,]), curve(TypeII.h(x, a=a, h=h), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIa") {
			with(as.list(coefs[4,]), curve(TypeIIIa.h(x, a=a, h=h), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			mtext(type)
		} else if (type=="TypeIIIb") {
			with(as.list(coefs[5,]), curve(TypeIIIb.h(x, a=a, h=h, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			mtext(type)
		} else if (type=="all") {
			with(as.list(coefs[1,]), curve(Type0(x, a=a), add=TRUE, col=lcols[1], lty=ltys[1], lwd=lwd))
			with(as.list(coefs[2,]), curve(TypeI(x, a=a), add=TRUE, col=lcols[2], lty=ltys[2], lwd=lwd))
			with(as.list(coefs[3,]), curve(TypeII.h(x, a=a, h=h), add=TRUE, col=lcols[3], lty=ltys[3], lwd=lwd))
			with(as.list(coefs[4,]), curve(TypeIIIa.h(x, a=a, h=h), add=TRUE, col=lcols[4], lty=ltys[4], lwd=lwd))
			with(as.list(coefs[5,]), curve(TypeIIIb.h(x, a=a, h=h, theta=theta), add=TRUE, col=lcols[5], lty=ltys[5], lwd=lwd))
			legend("bottomright",c("Type 0", "Type I","Type II","Type IIIa", "Type IIIb"), col=lcols, lty=ltys, lwd=lwd, bty="n")
			mtext("All fitted types")
		}
	}
}
