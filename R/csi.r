# Combined standardized indices

# Citations:
# de la Mare, W.K., and Constable, A.J. 2000. Utilising data from ecosystem monitoring for managing fisheries: 
# development of statistical summaries of indices arising from the CCAMLR ecosystem monitoring program. CCAMLR Sci. 7: 101-117.
# Boyd, I.L., and Murray, A.W.A. 2001. Monitoring a marine ecosystem using responses of upper trophic level predators. J. Anim. Ecol. 70: 747-760.
# Reid, K., Croxall, J.P., Briggs, D.R., and Murphy, E.J. 2005. 
# Antarctic ecosystem monitoring: quantifying the response of ecosystem indicators to variability in Antarctic krill. ICES J. Mar. Sci. 62: 366-373.
# See also chapter 11 by Croxall in
# Boyd, I., Wanless, S., and Camphuysen, C.J. 2006. Top Predators in Marine Ecosystems.

csi <- function(x) {
	# From Box 11.3 in Top Predators in Marine Ecosystems, eds. I.L.Boyd,S. Wanless and C. J. Camphuysen.
	# Croxall: Monitoring predator-prey interactions using multiple predator species: the South Georgia experience
	# Missing values? How to deal with them?
	
	x <- scale(x)
	a <- ifelse(!is.na(x), 1, 0)
	
	I <- sapply(1:nrow(x), function(i) t(a[i,]) %*% x[i,])
	# Edit: check that matrix is positive and semi-definite (some eigenvalues are negative)
	S <- cov(x)
	V <- diag(a %*% S %*% t(a))
	csi <- I / sqrt(V)
	list(csi=csi, I=I, V=V, S=S, eig=eigen(S), x=x)
}
