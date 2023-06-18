source("W:/PROJEKT/R/distributions/distRibutions.r")
source("DistributionsFunctions.r")

# Test the custom functions ecdf2 and ecdfn. Compare with the base function ecdf.

# Random sample
x <- rnorm(100,0,2)

# Built-in function ecdf(x)
ecdf.x <- ecdf(x)
dev.new(); plot(ecdf.x)

# Custom function ecdf2
test1 <- ecdf2(x)

# Matrix with equal sample sizes
# Same mean & n, vary sd
x1 <- rnorm(100,0,2)
x2 <- rnorm(100,0,4)
x3 <- rnorm(100,0,6)
x.mat <- cbind(x1,x2,x3)

test2 <- ecdfn(x.mat)

# List with unequal sample sizes
# Vary mean & n, same sd
x1 <- rnorm(100,0,2)
x2 <- rnorm(200,1,2)
x3 <- rnorm(500,2,2)
x.ls <- list(x1,x2,x3)

test3 <- ecdfn(x.ls)
