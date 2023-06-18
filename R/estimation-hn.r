# DISTANCE ESTIMATION
# Half-normal detection function (key without adjustment series)
# This closed-form estimator can only be used for the half-normal key,
# otherwise numerical solutions to M-L methods are necessary...
# Not exactly equal to DISTANCE-results (but very close)

# Only point estimation so far,
# no interval estimation

# Required input: an object with fields (columns) PerpDistance, ClusterSize
# Source: section 3.3.4 Buckland et al. 2001

dist.hn <- function(dat, w = FALSE) {
  # dat = input object with the data
  
  rm.ind <- which(dat$PerpDistance<0 | is.na(dat$Object) | is.na(dat$PerpDistance) | is.na(dat$ClusterSize) | dat$ClusterSize<0)
  
  x <- dat$PerpDistance[-rm.ind]
  cs <- dat$ClusterSize[-rm.ind]
  
  if (is.numeric(w)) w <- w
  if (w==FALSE) w <- max(dat$PerpDistance, na.rm=T)
  
  g <- function(x,sigma) exp(-(x^2)/(2*sigma^2))
  
  sumdist <- sum(x^2)
  nobs <- length(x)
  sigma <- sqrt(sumdist/nobs)
  
  gx <- g(x = x, sigma = sigma)
  # mu <- integrate(g, sigma=sigma, lower=0, upper=w)$value # gx dx
  mu <- integrate(g,
                  sigma = sigma,
                  lower = 0,
                  upper = Inf)$value # gx dx
  # mu <- sqrt((pi*sigma^2)/2)
  fx <- gx / mu
  
  #f0 <- 1 / mu
  f0 <- sqrt(2/(pi*sigma^2))
  esw <- mu
  
  xv <- seq(from = 0, to = w, by = 0.1)
  gxv <- g(x = xv, sigma = sigma)
  fxv <- gxv / mu
  
  par(mfrow=c(1,2))
  plot(xv, gxv, font.lab=2, las=1, type="l", ylab="g(x)", main = "Detection function")
  plot(xv, fxv, font.lab=2, las=1, type="l", ylab="f(x)", main = "Probability density function")
  abline(h=f0, lty=2, col=2)
  par(mfrow=c(1,1))
  
  # Effort
  l <- sum(as.vector(tapply(dat$LineLength,dat$LineLabel,sum)/tapply(dat$LineLength,dat$LineLabel,length)))
  samples <- length(unique(dat$LineLabel))
  
  a <- 2*w*l
  p <- mu / w
  
  # Density of clusters
  ds <- 1000*(nobs*f0)/(2*l)
  
  # es (expected cluster size)
  avg.s <- mean(cs)
  # Calculate size biased estimate of expected cluster size
  # see section 3.5.4 in Distance book
  xv <- hn(x = x, sigma = sigma)
  z <- log(cs)
  fm <- lm(z ~ xv)
  es <- exp(sum(coef(fm)) + var(z)/2)
  
  # Density
  d <- 1000*(nobs*f0*es) / (2*l)
  
  # Short output (compare with DISTANCE)
  out <- list(
    rm = rm.ind,
    Effort = l,
    Samples = samples,
    Width = w,
    Observations = nobs,
    sigma = sigma,
    f0 = f0,
    ESW = esw,
    p = p,
    area = a,
    DS = ds,
    mean.S = avg.s,
    E.S = es,
    D = d)
  
  print(out)
}
