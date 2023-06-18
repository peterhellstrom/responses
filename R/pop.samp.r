# Custom functions for standard error & coefficient of variation
se <- function(x) sd(x) / sqrt(length(x))
cv <- function(x) sd(x) / mean(x)
ci.l <- function(x) {
	mean(x) - qnorm(0.975)*se(x)
}
ci.u <- function(x) {
	mean(x) + qnorm(0.975)*se(x)
}

# ToDO: Handling of NA's in tapply.formula?
pop.samp <- function (fo, df, rownames=TRUE, ci=FALSE) {
	
	df <- droplevels(df)
	z <- tapply.formula(fo, df, length)
	dims <- dim(z)
	out.names <- dimnames(z)
	
	if (ci == TRUE) {
		funcs <- c("sum", "length", "mean", "median", "sd", "se", "cv", "min", "max", "ci.l", "ci.u")
		out <- sapply(1:length(funcs), function (i) do.call(tapply.formula, list(fo=fo, df=df, func=get(funcs[i]))))
		colnames(out) <- c("Sum", "n", "Mean", "Median", "SD", "SE", "CV", "min", "max", "CIL", "CIU")
	} else if (ci == FALSE) {
		funcs <- c("sum", "length", "mean", "median", "sd", "se", "cv", "min", "max")
		out <- sapply(1:length(funcs), function (i) do.call(tapply.formula, list(fo=fo, df=df, func=get(funcs[i]))))
		colnames(out) <- c("Sum", "n", "Mean", "Median", "SD", "SE", "CV", "min", "max")
	}

	if (length(dims) > 1) {
		g <- expand.grid(out.names)
		g.row <- apply(g[,1:ncol(g)], 1, paste, collapse = ".")
		out <- data.frame(out, g)
		if (rownames) rownames(out) <- g.row
		out
	} else {
		out <- as.data.frame(out)
		out <- data.frame(out, out.names)
		if (rownames == FALSE) rownames(out) <- NULL
		out
	}
}
