# Common detection functions
# hn = half-normal
# hr = hazard rate
# u = uniform
# Should add negative exponential for completeness

hn <- function(x,sigma) exp(-(x^2)/(2*sigma^2))
hr <- function(x,sigma,b) 1 - exp(-(x/sigma)^-b)
