# A simple grid search function
# This is one way to find start-values for estimation of functional response parameters.
# NOTE: the nls2 package has more advanced grid search options for start values!

#' @export
fr.start.grid <- function(x, y, type=c("TypeII","TypeIIIa","TypeIIIb"), length.out=10, extra=50) {

	type <- match.arg(type)

	# Generate grid
	# For
	# a (asymptote): search in interval mean(y) to max(y) + extra
	# b (half-saturation constant): 0 to mean(x)
	# theta: 0.5 to 5

	if (type == "TypeII") {
		grid <- expand.grid(a = seq(mean(y), max(y) + extra, length.out=length.out), b = seq(0, mean(x), length.out=length.out))
		ss <- function(p) sum((y - TypeII(x, p[1], p[2]))^2)
		idx  <- which.min(apply(grid, 1, ss))
		startval <- grid[idx,]
		fm <- nls(y ~ TypeII(x, a, b), start = startval)
	} else if (type == "TypeIIIa") {
		grid <- expand.grid(a = seq(mean(y), max(y) + extra, length.out=length.out), b = seq(0, mean(x), length.out=length.out))
		ss <- function(p) sum((y - TypeIIIa(x, p[1], p[2]))^2)
		idx  <- which.min(apply(grid, 1, ss))
		startval <- grid[idx,]
		fm <- nls(y ~ TypeIIIa(x, a, b), start = startval)
	} else if (type == "TypeIIIb") {
		grid <- expand.grid(a = seq(mean(y), max(y) + extra, length.out=length.out), b = seq(0, mean(x), length.out=length.out), theta = seq(0.5, 5, length.out=length.out))
		ss <- function(p) sum((y - TypeIIIb(x, p[1], p[2], p[3]))^2)
		idx  <- which.min(apply(grid, 1, ss))
		startval <- grid[idx,]
		fm <- nls(y ~ TypeIIIb(x, a, b, theta), start = startval)
	} else stop("Model not included!")

	dat <- list('start'=startval, 'nls'=fm, 'nlsSum'=summary(fm))
}
