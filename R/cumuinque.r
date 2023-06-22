# Calculate the cumulative number of unique values in a vector
#' @export
cumunique <- function(x) cumsum(ifelse(!duplicated(x),1,0))
