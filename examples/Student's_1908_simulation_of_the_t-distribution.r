# http://facultypages.morris.umn.edu/~jongmink/Stat2611/s1.pdf

n <- 3000
n_groups <- 750

mu <- 10
sigma <- 1

x <- rnorm(n, mu, sigma)
g <- rep(1:n_groups, each = n / n_groups)

# Calculate t-statistic

tst <- (tapply(x, g, mean) - mean(x)) / (tapply(x, g, sd) / sqrt(tapply(x, g, length)))

plot(density(tst), main = "Simulation of Student's t")
lines(density(scale(x)), col = 3)

curve(dnorm(x, 0, sigma), n = 1001, col = 2, add = TRUE)
curve(dt(x, df = (n / n_groups) - 1), col = 4, add = TRUE)

legend("topright", c("Normal", "t"), lty = c(1, 1), col = c(2, 4), bty = "n")
