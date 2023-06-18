# Inverse logit
expit <- function(x) {
  exp(x) / (1+exp(x))
}

# logit
logIt <- function(x) {
  1/(1+1/(exp(x)))
}
