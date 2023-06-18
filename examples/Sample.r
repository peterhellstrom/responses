##################################################
# Create random sample, unweighted

xv <- 1:5
n <- 100000
.x <- sample(xv, replace=TRUE, size=n)

table(.x) / n # should be close to 1/length(xv) for all categories
1/length(xv)

hist(.x, freq=FALSE)

# Leave out argument n:
# Return a sample of length xv
.x <- sample(xv, replace=TRUE)

##################################################
# Weighted random sample, use argument 'prob' as weights

xv <- 1:5
n <- 100000
# Probabilities doesn't have to sum to one (but have to be non-negative and not all zero)
ps <- c(300,400,500,700,600)
# R recalculates the prob vector internally, probabilites should in this case be:
ps/sum(ps)

length(xv) == length(ps) # Must evaluate to TRUE

.x <- sample(xv, replace=TRUE, size=n, prob=ps)

# Compare empirical with theoretical sample probabilities
table(.x) / n
ps/sum(ps)

hist(.x, freq=FALSE)

##################################################

# library(sampling) has more advanced functions.
