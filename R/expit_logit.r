# Inverse logit
#' @export
expit <- function(x) {
  exp(x) / (1+exp(x))
}

# logit
#' @export
logIt <- function(x) {
  1/(1+1/(exp(x)))
}
