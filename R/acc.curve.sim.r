# Restructure code like this instead.
# sim.gen.data <- function(genotypes, samples, resamplings) ()
# sim.run.data <- function(input.data) ()
# sim.summary <- function(run.data) ()
# sim.plot <- function(sim.summary) ()

#' @export
acc.curve.sim <- function(genotypes, samples, resamplings) {

	# Inputs: genotypes (list with true population sizes), samples, resamplings

	pop.vec <- sapply(1:length(genotypes), function(i) max(genotypes[[i]]) )

	# Vectors to hold output (g is for global)
	avg.output <- vector("list",length=length(samples))
	full.output <- vector("list",length=length(samples))
	g.avg.output <- vector("list",length(genotypes))
	g.full.output <- vector("list",length(genotypes))
	avg.prop.found <- matrix(NA,ncol=length(genotypes),nrow=length(samples))
	colnames(avg.prop.found) <- pop.vec
	rownames(avg.prop.found) <- samples
	full.prop.found <- matrix(NA,ncol=resamplings,nrow=length(samples))
	g.full.prop.found <- vector("list",length(genotypes))

	# Outer loop - loops over length of true pop size (input genotypes)
	# Inner loop - sample size

	for (runs in 1:length(genotypes)) {

		for (res in 1:length(samples)) {

		# Create data sets
		x <- 1:samples[res]
		y <- matrix(NA,ncol=resamplings,nrow=samples[res])
			for (j in 1:resamplings) y[,j]<- sample(genotypes[[runs]],samples[res],replace=TRUE)

		z <- matrix(NA,ncol=resamplings,nrow=samples[res])
			for (k in 1:resamplings) z[,k] <- rep(1:length(unique(y[,k])),as.vector(table(y[,k])) )
			z.mat <- matrix(rep(z,length(resamplings)),nrow=samples[res],ncol=resamplings,byrow=FALSE)
			shuffle <- apply(z.mat,2,sample,replace=FALSE)

		data.set <- matrix(NA,nrow=dim(shuffle)[1],ncol=dim(shuffle)[2])
			for (l in 1:ncol(shuffle)) {
				for (m in 1:nrow(shuffle)) data.set[m,l] <- length(unique(shuffle[1:m,l]) ) }


		# Extract proportion of found genotypes (averaged over resamplings)
		a <- max(genotypes[[runs]])
		b <- nrow(data.set)
		avg.prop.found[res,runs] <- ("Mean prop found"=as.vector(summary(data.set[b,])[4])/a )

		# Extract proportion found for each individual run
		full.prop.found[res,] <- data.set[b,]/a

		# Non-linear curve fitting
		kohn <- lapply(1:resamplings, function(i) try(summary(nls(y~SSmicmen(x,a,b),data=list(x=x,y=data.set[,i]))) ) )

		kohn <- kohn[sapply(kohn, function(x) !inherits(x, "try-error"))]
		kohn.l <- length(kohn)
		kohn.est <- t(sapply(1:kohn.l, function(i) kohn[[i]]$parameters[1,]))

		eggert <- lapply(1:resamplings, function(i) try(summary(nls(y~SSasympOrig(x,a,b),data=list(x=x,y=data.set[,i]))) ) )
		eggert <- eggert[sapply(eggert, function(i) !inherits(i, "try-error"))]
		eggert.l <- length(eggert)
		eggert.est <- t(sapply(1:eggert.l, function(i) eggert[[i]]$coefficients[1,]))

		chessel <- lapply(1:resamplings, function(i) try(summary(nls(y~Chessel(x,a),start=list(a=median(eggert.est[,1])),data=list(x=x,y=data.set[,i])) ) ) )
		chessel <- chessel[sapply(chessel, function(i) !inherits(i, "try-error"))]
		chessel.l <- length(chessel)
		chessel.est <- t(sapply(1:chessel.l, function(i) chessel[[i]]$coefficients[1,]))

	# Create output with average values (stored in avg.output)
	quants <- function (x) as.vector(quantile(x,c(0.025,0.975) ) ) # Get 95% quantile interval
	# Gather quantile intervals in a matrix
	qi <- matrix(c(quants(kohn.est[,1]),quants(eggert.est[,1]),quants(chessel.est[,1])),ncol=3,nrow=2)
	# Gather number of succesful runs
	coef.l <- as.vector(c(length(kohn.est[,1]),length(eggert.est[,1]),length(chessel.est[,1])))
	# Combine parameter estimates & intervals in one output
	avg.out <- rbind(cbind(as.vector(summary(kohn.est[,1])), as.vector(summary(eggert.est[,1])), as.vector(summary(chessel.est[,1])) ),qi,coef.l)
	# Set row- & colnames
	rownames(avg.out) <- c("Min","1st Qu.","Median","Mean","3rd Qu.","Max.","2.5% Qu.","97.5% Qu.","N")
	colnames(avg.out) <- c("kohn","eggert","chessel")
	# Write avg.output to list-object
	avg.output[[res]] <- avg.out
	full.output[[res]] <- list(kohn.est=kohn.est,eggert.est=eggert.est,chessel.est=chessel.est)
	}
	g.avg.output[[runs]] <- avg.output
	g.full.output[[runs]] <- full.output
	g.full.prop.found[[runs]] <- full.prop.found
	}
	list(average=g.avg.output,full=g.full.output,AvgPropFound=avg.prop.found,FullPropFound=g.full.prop.found)
}
