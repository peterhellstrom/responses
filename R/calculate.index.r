# REMOVAL TRAPPING INDEX
# Index for relative abundance
# Various methods compiled from the literature

# Terminology (compared to the ML functions):
# a2 = captures
# T = effort (in ml.ce)

# Input values can be single values, or a matrix/list(?) containing data for many
# trapping sessions.

# Requires four inputs, nontargets (a1), targets (a2), traps (T) & duration (D)
# T = number of traps
# D = duration of study, number of days (total uncorrected trap effort = T * D)
# a0 = number of traps that were not sprung
# a1 = sum of the number of sprung, but empty, traps and traps that caught a non-target species
# a2 = number of traps that caught an individual of the target species

# Index from Fitzgerald, Efford & Karl
# New Zealand Journal of Zoology, 31:167-184 2004
# If a1 = 0, the point estimate of Fitzgerald's index equals Caughley's index (see below)
# Estimation of SE assumes that captures per trap follow a Poisson distribution
# 2) Compare Fitzgerald's estimate to Caughley's (from his classic book, 1977)
# Caughley did not account for competing species, only the target species
# If a1 in Fitzgerald's index > 0, then Caughley's estimator is biased low
# NOTE: Fitzgerald's index is really equation 4 in Linn & Downton (1975)
# 3) Nelson & Clark's (1973) index (Journal of Mammalogy, 1973)
# This index takes number of checks/total trap time into account
# The formula below is a modification by Peter Hellstr?m
# New, additional parameter:
# cpt = checks per trap (number of times each trap was checked)
# 4) Uncorrected density index
# Standard approach in many studies

DI <- function(targets, nontargets=0, traps, duration, cpt=duration, rowNames=NULL){
	
		a1 <- nontargets
		a2 <- targets
		effort <- traps*duration
		correction <- duration*(a1 + a2)/(2*cpt)
		corrected.effort <- effort - correction
		a0 <- effort - a1 - a2
	
		index <- data.frame(
					fitzgerald = -100*(a2/(a1+a2))*log(a0/effort),
					caughley = -100*log(1-(a2/effort)),
					nelson = 100*a2/corrected.effort,
					uncorrected = 100*a2/effort
				)
	
		# Fitzgerald
		se.index <- data.frame(
						fitzgerald = 100*sqrt(((a1*a2*(log(a0/effort)^2))/(effort-a0)^3)+(a2^2/(a0*effort*(effort-a0))))
					)
		
		input <- data.frame(
					targets = targets,
					nontargets = nontargets,
					traps = traps,
					duration = duration,
					cpt = cpt,
					effort = effort,
					correction = correction,
					corrected.effort = corrected.effort)

		if (is.character(rowNames)==TRUE | is.numeric(rowNames)) {
			if (length(rowNames) != nrow(index)) stop("Rownames does not match")
			rownames(input) <- rownames(index) <- rownames(se.index) <- rownames
		}
		
		out <- list(
			input = input,
			index = index,
			se.index = se.index
			)
		
		out
	}
