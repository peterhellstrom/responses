acc.sim.bias.ratio.plot <- function(genotypes, samples, object, model) {

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )
	.x <- sapply(1:length(genotypes), function(i) samples/pop.vec[i])
 
	model.lab <- model

	if (model=="kohn") model <- 1
	else if (model=="eggert") model <- 2
	else if (model=="chessel") model <- 3
	else stop ("Wrong model chosen!")

	xlims <- max(.x)
	ylims.lo <- min(unlist(object[[model]]$bias))
	ylims.up <- max(unlist(object[[model]]$bias))

	dev.new(width=6,height=5)
	par(mar=c(5,4,3,1))
	plot(c(0,xlims),c(ylims.lo,ylims.up),type="n",xlab="Sample size / True population size",ylab="Relative bias",font.lab=2,las=1,main=model.lab)

	# Plot point estimates
	xv <- as.vector(.x)
	yv <- as.vector(object[[model]]$bias)
	cols <- rep(1:length(genotypes), each=length(samples))
	points(xv,yv,pch=16,col=cols,cex=1.3)

	# Add legend
	legend("topright",legend=pop.vec, col=c(1:length(genotypes)),pch=16,cex=1.3 )
}
