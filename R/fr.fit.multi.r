# Multi-prey
# This code only estimates fr to one prey species,
# update to multi-response function (LLmulti.disc 2 or similar)

fr.fit.multi <- function(data) {

x1 <- data$x1
x2 <- data$x2
y1 <- data$y1
y2 <- data$y2

nlc <- nls.control(maxiter = 50000, tol=0.000001, minFactor=0.000001)

  inds <- which(x1<=0 | x2<=0 | y1<=0 | y2<=0)
  
  if (length(inds) > 0) {
  
  y1.new <- y1[-inds]
  y2.new <- y2[-inds]
  x1.new <- x1[-inds]
  x2.new <- x2[-inds]  

  mm.start.mle1 <- as.list(getInitial(y1.new ~ SSmicmen(x1.new, a, b),
              data=data.frame(x1.new,y1.new), control=nlc))
  start.mle1 <- list(a1=as.numeric((mm.start.mle1$a/mm.start.mle1$b)), 
              h1=as.numeric(1/mm.start.mle1$a))
                            
  mm.start.mle2 <- as.list(getInitial(y2.new ~ SSmicmen(x2.new, a, b),
              data=data.frame(x2.new,y2.new), control=nlc))
  start.mle2 <- list(a2=as.numeric((mm.start.mle2$a/mm.start.mle2$b)), 
              h2=as.numeric(1/mm.start.mle2$a))
              
  fm.multi <- mle2(minuslogl=LLmulti.disc2bd, start=
              list(a1=start.mle1$a1, a2=start.mle2$a2, h1=start.mle1$h1, h2=start.mle2$h2), 
              data=list(x1=x1,x2=x2,y1=y1,y2=y2), control=list(maxit=50000))
  }
  
  if (length(inds) == 0) {
  
  mm.start.mle1 <- as.list(getInitial(y1 ~ SSmicmen(x1, a, b),
              data=data.frame(x1,y1), control=nlc))
  start.mle1 <- list(a1=as.numeric((mm.start.mle1$a/mm.start.mle1$b)), 
              h1=as.numeric(1/mm.start.mle1$a))
                            
  mm.start.mle2 <- as.list(getInitial(y2 ~ SSmicmen(x2, a, b),
              data=data.frame(x2,y2), control=nlc))
  start.mle2 <- list(a2=as.numeric((mm.start.mle2$a/mm.start.mle2$b)), 
              h2=as.numeric(1/mm.start.mle2$a))
              
  fm.multi <- mle2(minuslogl=LLmulti.disc2bd, start=
              list(a1=start.mle1$a1, a2=start.mle2$a2, h1=start.mle1$h1, h2=start.mle2$h2), 
              data=list(x1=x1,x2=x2,y1=y1,y2=y2), control=list(maxit=5000))
  }
  
coefs <- rbind(
coef(fm.multi)
)

# colnames(coefs) <- c("a1","a2","h1","h2")
rownames(coefs) <- c("MultiDisc")

# Calculate AICc and AICc-weights:
# Function that calculates corrected AIC:
AICc <- function(object) {
k <- length(coef(object)) + 1
n <- length(x1)
AICc <- as.numeric(-2*logLik(object) + 2*k) + (2*k*(k+1)/(n-k-1))
out <- c(logLik(object),k,AICc)
names(out) <- c('LogL','K','AICc')
out
}

AICc.table <- c(
as.vector(AICc(fm.multi)[3])
)
names(AICc.table) <- c("MultiDisc")

deltai <- AICc.table - min(AICc.table)
rel.like <- exp(-deltai/2)
wi <- rel.like / sum(rel.like)
names(wi) <- c("MultiDisc")

fit.objects <- list(fm.multi=fm.multi)

out <- list(x1=x1, x2=x2, y=y1, start.mle1=start.mle1, start.mle2=start.mle2, 
coefs=coefs, AICc=AICc.table, AkaikeWeights=wi, fit.objects=fit.objects)
out

}
