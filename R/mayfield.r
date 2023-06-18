# Mayfield - estimator
# Flint et al 1995
#######################################################################

mayfield <- function(dat, method) {

	n.broods <- dim(dat)[1]
	cols <- ncol(dat)

	if (is.null(ncol(dat)) == TRUE) stop("Data frame must have at least two columns!")
	n.obs <- apply(dat, 1, function(x) length(which(!is.na(x))))
	n.intervals <- n.obs - 1

	# Calculate delta brood survival
	# Find all observations on a row
	inds <- lapply(1:n.broods, function(i) which(!is.na(dat[i,])))
	dat.obs <- lapply(1:n.broods, function(i) dat[i,inds[[i]]])

	# Calculate brood size at first observation - brood size at final observation
	delta.bs <- sapply(1:n.broods, function(i) as.numeric(dat.obs[[i]][1] - dat.obs[[i]][length(dat.obs[[i]])]))

	time <- unlist(strsplit(colnames(dat),"X"))
	time <- as.numeric(time[time!=""])
	time.obs <- lapply(1:n.broods, function(i) time[inds[[i]]])

	# If observations are made at equally spaced intervals:
	L <- sapply(1:length(time.obs), function(i) unique(diff(time.obs[[i]])))
	n.L <- length(unique(n.obs))

	###
	if (method == "exposure-days") {

		if (n.L == 1) {
			cat(paste("Broods were observed at equally spaced intervals, L = ", unique(L)), "\n")
			if (cols == 2) sr <- dat[,-1] / dat[,-ncol(dat)]
			if (cols > 2) sr <- rowSums(dat[,-1]) / rowSums(dat[,-ncol(dat)])
			dsr <- sr^(1/L)
			exposure <- delta.bs / (1-dsr)
			exposure2 <- rowSums(L*(dat[-1] + dat[-ncol(dat)])*0.5)
			# exposure not defined if delta.bs == 0, and dsr == 1. Returns NaN
			# Replace with values from exposure2
			nan.inds <- which(is.nan(exposure))
			exposure[nan.inds] <- exposure2[nan.inds]
		}

		if (n.L > 1) {
			cat("Broods were not observed at equally spaced intervals", "\n")
			# Calculate exposure for each brood
			# Step 1: Calculate sum of two consecutive counts / broods
			step1 <- lapply(1:n.broods, function(i) dat.obs[[i]][-length(dat.obs[[i]])] + dat.obs[[i]][-1])
			# Step 2: Calculate length of each interval
			step2 <- lapply(1:n.broods, function(i) diff(time.obs[[i]]))
			# Step 3: Calculate exposure (assuming changes in brood size occur at the midpoint of the observation interval)
			step3 <- lapply(1:n.broods, function(i) step1[[i]] * step2[[i]] * 0.5)
			# Final step, sum all exposure days
			exposure <- sapply(1:n.broods, function(i) sum(step3[[i]]))
			# Daily Survival Rate
			dsr <- 1 - (delta.bs / exposure)
		}

	} else if (method == "bart-robson") {

		pm.mayfield <- function(alive,dead,intrvl) 1 - (sum(dead) / sum((intrvl*(alive + dead*0.5))))
		dsr <- numeric(length(dat.obs))
		exposure <- numeric(length(dat.obs))

		for (i in 1:length(dat.obs)) {
			temp <- as.numeric(dat.obs[[i]])
			dead <- temp[-length(temp)] - temp[-1]
			alive <- temp[-length(temp)] - dead
			intrvl <- diff(time.obs[[i]])

			# Original estimate
			p.hat <- pm.mayfield(alive,dead,intrvl)
			if (p.hat == 1 & sum(dead)==0) {
				dsr[i] <- p.hat
			} else {
				n.iter <- 10
				p.iter <- numeric(n.iter)

				for (j in 1:n.iter) {
					f.pm <- sum((intrvl/p.hat) * (alive - (dead*p.hat^intrvl / (1 - p.hat^intrvl))))
					f.prim.pm <- sum((intrvl/p.hat^2) * (alive + ((dead*p.hat^intrvl)*(intrvl-1+p.hat^intrvl)) / (1-p.hat^intrvl)^2))
					p.hat <- p.hat + f.pm/f.prim.pm
					p.iter[j] <- p.hat
				} # end of inner loop
	
				dsr[i] <- p.iter[n.iter]
				# Problem here: exposure is not defined if dsr == 1, and prints as zero.
				exposure[i] <- delta.bs[i] / (1 - dsr[i])
			} # end of else

		} # end of outer loop

		# Calculate exposure2 for each brood
		# Step 1: Calculate sum of two consecutive counts / broods
		step1 <- lapply(1:n.broods, function(i) dat.obs[[i]][-length(dat.obs[[i]])] + dat.obs[[i]][-1])
		# Step 2: Calculate length of each interval
		step2 <- lapply(1:n.broods, function(i) diff(time.obs[[i]]))
		# Step 3: Calculate exposure (assuming changes in brood size occur at the midpoint of the observation interval)
		step3 <- lapply(1:n.broods, function(i) step1[[i]] * step2[[i]] * 0.5)
		# Final step, sum all exposure days
		exposure2 <- sapply(1:n.broods, function(i) sum(step3[[i]]))

		nan.inds <- which(exposure==0)
		exposure[nan.inds] <- exposure2[nan.inds]

	} # end of method bart-robson

	###
	# Calculate standard error
	se.dsr <- function(exposure, dsr) {
		M <- n.broods
		exposure.bar <- sum(exposure/M)
		dsr.bar <- 1 - (sum(delta.bs) / sum(exposure))
		num <- sum((exposure^2)*((dsr - dsr.bar)^2))
		denom <- ((exposure.bar^2)*(M - 1)*M)
		sqrt(num / denom)
	}

	dsr.pt <- 1 - (sum(delta.bs) / sum(exposure))
	dsr.int <- se.dsr(exposure,dsr)

	###
	# Collect output

	out <- list(
		# List with brood level output
		Brood = data.frame(n.obs, delta.bs, exposure, dsr),
		# Vector with population level output
		PointEstimate = c(n = sum(n.obs), delta.bs = sum(delta.bs), exposure = sum(exposure), dsr = dsr.pt),
		# Vector with standard error and confidence interval
		IntervalEstimate = c(se = dsr.int, cil = dsr.pt - 2*dsr.int, ciu = dsr.pt + 2*dsr.int)
	) # End of output list

	# Print output
	out

}
