# Calculate the cumulative number of unique values in a vector
cumunique <- function(x) cumsum(ifelse(!duplicated(x),1,0))
