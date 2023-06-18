# Plot of average values
# (Estimated pop size against sample size)

acc.sim.avg.plot <- function(genotypes, samples, object, model, type) {

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )

	model.lab <- model

	if (model=="kohn") model <- 1
	else if (model=="eggert") model <- 2
	else if (model=="chessel") model <- 3
	else stop ("Wrong model chosen!")

	if (type=="qi") {
		xlims <- max(object[[model]]$samples) 
		ylims <- max(unlist(object[[model]]$qi))

	dev.new(width=6,height=5)
	par(mar=c(5,4,3,1))
	plot(c(0,xlims),c(0,ylims),type="n",xlab="Sample size",ylab="Estimated population size",font.lab=2,las=1,main=model.lab)

	# Plot point estimates
	xv <- rep(object[[model]]$samples,times=length(genotypes) )
	yv <- as.vector(object[[model]]$avg.est)
	cols <- rep(1:length(genotypes), each=length(samples))
	points(xv,yv,pch=16,col=cols,cex=1.3)

	# Plot quantile intervals
	for (i in 1:length(genotypes)) {
		xvec <- samples
		low.arr <- object[[model]]$qi[[i]][,1]
		up.arr <- object[[model]]$qi[[i]][,2]
		arrows(xvec,low.arr,xvec,up.arr,angle=90,code=3,length=0.05) }

	# Plot true population size as horizontal line
	pop.vec <- colnames(object[[model]]$avg.est)
	for (i in 1:length(genotypes)) {
		abline(h=pop.vec[i],lty=2,col=c(0+i))
	}} 

	if (type=="point") {
		xlims <- max(object[[model]]$samples) 
		ylims <- max(unlist(object[[model]]$avg.est))

		dev.new(width=6,height=5)
		par(mar=c(5,4,3,1))
		plot(c(0,xlims),c(0,ylims),type="n",xlab="Sample size",ylab="Estimated population size",font.lab=2,las=1,main=model.lab)

		# Plot point estimates
		xv <- rep(object[[model]]$samples,times=length(genotypes) )
		yv <- as.vector(object[[model]]$avg.est)
		cols <- rep(1:length(genotypes), each=length(samples))
		points(xv,yv,pch=16,col=cols,cex=1.3)

		# Plot true population size as horizontal line
		pop.vec <- colnames(object[[model]]$avg.est)
		for (i in 1:length(genotypes)) {
		abline(h=pop.vec[i],lty=2,col=c(0+i))
		}}

	# Add legend
	legend("topright",legend=pop.vec, col=c(1:length(genotypes)),pch=16,cex=1.3 ) 
}
