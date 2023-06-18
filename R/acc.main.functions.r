# This R-code uses three estimation functions that have been used in the literature

Kohn <- function(x,a,b) a*x/(b+x) # Michaelis-Menten
Eggert <- function(x,a,b) a*(1-exp(b*x)) # Two-parameter exponential
Chessel <- function(x,a) a-(a*(1-(1/a))^x)

Kohn.gr <- deriv(~ a*x / (b + x), c("a","b"), function(a,b,x) {})
Eggert.gr <- deriv(~ a*(1-exp(b*x)), c("a","b"), function(a,b,x) {})
Chessel.gr <- deriv(~ a-(a*(1-(1/a))^x), c("a"), function(a,x) {})

# Note that both the Michaelis-Menten equation and the two-parameter exponential
# equations are available in R as SSmicmen resp. SSasympOrig.
