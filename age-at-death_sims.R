library(data.table)
library(plyr)
library(devtools)
# [LOAD archSeries]

# Define breaks (just for axis labels)
breaks <- c(0,2,6,12,24,36,48,72,96,120)

# Define models in terms of survivorship [NEED TO ADD MORE]
security = c(0.904, 0.904, 0.645, 0.38, 0.25, 0.239, 0.171, 0.118, 0)
milk = c(0.47, 0.42, 0.39, 0.35, 0.28, 0.23, 0.19, 0.1, 0)
wool = c(0.85, 0.75, 0.65, 0.63, 0.57, 0.5, 0.43, 0.2, 0)

# Convert to mortality
security <- -diff(c(1, security))
milk <- -diff(c(1, milk))
wool <- -diff(c(1, wool))

# Define function to simulate KS test results
sim.sig <- function(n, models,  probs1, probs2, reps=1000, test="ks") {
  
  #Set up bin labels
  labels <- numeric(0)
  for(i in 1:(length(breaks) - 1)) {
    labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
  }
  
  #Perform simulation
  rep.no <- rep(1:reps, each = n)
  sim <- data.table(rep.no=rep(1:reps, each=n), id=rep(1:n, reps))
  sim[, a := sample(1:length(labels), size=nrow(sim), replace=TRUE, prob=probs1)]
  sim[, b := sample(1:length(labels), size=nrow(sim), replace=TRUE, prob=probs2)]
  
  #Calculate totals per bin
  totals1 <- sim[, j=list(count.a=length(id)), by=list(rep.no, a)]
  totals2 <- sim[, j=list(count.b=length(id)), by=list(rep.no, b)]
  setnames(totals1, "a", "b")
  totals <- merge(totals1, totals2, by=c("rep.no", "b"), all=TRUE)
  
  #Tidy up results
  frame <- data.table(rep.no=rep(1:reps, each=length(labels)), b=rep(1:length(labels), reps), bin=rep(labels, reps))
  totals <- merge(frame, totals, by=c("rep.no", "b"), all=TRUE)
  totals[is.na(totals)] <- 0
    
  #Convert to cumulative
  cum.totals <- totals[, j=list(a=cumsum(count.a), b=cumsum(count.b)), by="rep.no"]
  cum.totals[, D := (a - b) / n]
  setnames(cum.totals, c("a", "b"), models)
  
  #Kolmogorov-Smirnov
  D.max <- cum.totals[, j=list(Dmax=max(abs(D))), by="rep.no"]
  cut.off <- 1.36 * sqrt(n*2/n^2)
  D.max[, sig := 0]
  D.max[Dmax > cut.off, sig := 1]
  x <- sum(D.max$sig)/reps
  
  #Chi-squared
  X2.fun <- function(totals) {
    x <- chisq.test(totals[,c("count.a", "count.b")])
    x[[3]]
  }
  X2 <- data.table(p=daply(.data = totals, .variables = "rep.no", .fun = X2.fun))
  X2 <- X2[!is.na(p),]
  X2[, sig := 0]
  X2[p < 0.05, sig := 1]
  y <- sum(X2$sig)/nrow(X2)
  
  #Return results
  c(x,y)
  
}

#Set up for simulation [here to allow testing within ks.sim - REMOVE LATER]
reps <- 1000
n <- 50
models <- c("wool", "security")
probs1 <- wool
probs2 <- wool

# [DEMONSTRATE CONCEPT AND CONFIDENCE INTERVALS USING poly.chron]

# Testing for false negatives
sizes <- seq(10, 250, 5)

false.negs <- 100 * (1 - sapply(sizes, sim.sig, models=c("wool", "security"), probs1=wool, probs2=security, reps=1000))
plot(sizes, false.negs[1,], main="Wool vs. Security", ylab = "False negative rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.negs[2,], ylab = "False negative rate (%)", xlab = "Sample size", col="darkblue", pch=19)
legend("topright", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.negs <- 100 * (1 - sapply(sizes, sim.sig, models=c("milk", "security"), probs1=milk, probs2=security, reps=1000))
plot(sizes, false.negs[1,], main="Milk vs. Security", ylab = "False negative rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.negs[2,], ylab = "False negative rate (%)", xlab = "Sample size", col="darkblue", pch=19)
legend("topright", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.negs <- 100 * (1 - sapply(sizes, sim.sig, models=c("milk", "wool"), probs1=milk, probs2=wool, reps=1000))
plot(sizes, false.negs[1,], main="Milk vs. Wool", ylab = "False negative rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.negs[2,], ylab = "False negative rate (%)", xlab = "Sample size", col="darkblue", pch=19)
legend("topright", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

# Testing for false positives
false.pos <- 100 * sapply(sizes, sim.sig, models=c("wool", "wool"), probs1=wool, probs2=wool, reps=10000)
plot(sizes, false.pos[1,], main="Wool vs. Wool", ylim=c(0,6.5), ylab = "False positive rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.pos[2,], ylab = "False positive rate (%)", xlab = "Sample size", col="darkblue", pch=19)
abline(h=5)
legend("topleft", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.pos <- 100 * sapply(sizes, sim.sig, models=c("milk", "milk"), probs1=milk, probs2=milk, reps=10000)
plot(sizes, false.pos[1,], main="Milk vs. Milk", ylim=c(0,6.5), ylab = "False positive rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.pos[2,], ylab = "False positive rate (%)", xlab = "Sample size", col="darkblue", pch=19)
abline(h=5)
legend("topleft", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.pos <- 100 * sapply(sizes, sim.sig, models=c("security", "security"), probs1=security, probs2=security, reps=10000)
plot(sizes, false.pos[1,], main="Security vs. Security", ylim=c(0,6.5), ylab = "False positive rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.pos[2,], ylab = "False positive rate (%)", xlab = "Sample size", col="darkblue", pch=19)
abline(h=5)
legend("topleft", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")


  

