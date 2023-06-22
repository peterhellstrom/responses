# After the simulation:
# Alternative power-analysis, check the number of times each model was selected
# as "best", using AICc as criterion:
# data is an object fitted with fr.fit or fr.fit.n
#' @export
fr.fit.out <- function (data){
	# A function that finds the position of the maximum value in a vector
	max.AICc <- function(x) which.max(x)
	# Calculate in which column the maximum value is found, and make a table:
	tab1 <- table(apply(data$AICc.weights,1,max.AICc))
	# Make a better table
	inds <- as.numeric(names(tab1))
	tab2 <- numeric(5)
	names(tab2) <- c("Type0","TypeI","TypeII","TypeIIIa","TypeIIIb")
	tab2[inds] <- tab1
	# Calculate "Power", percentage of simulations that each model was selected
	pow <- round(tab2/sum(tab2),2)
	# Output
	list(n=sum(tab2),table=tab2,power=pow)
}
