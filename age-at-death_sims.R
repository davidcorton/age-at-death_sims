library(data.table)

breaks <- c(0,2,6,12,24,36,48,72,96,120)

Security = c(0.904, 0.904, 0.645, 0.38, 0.25, 0.239, 0.171, 0.118, 0)
Milk = c(0.47, 0.42, 0.39, 0.35, 0.28, 0.23, 0.19, 0.1, 0)
Wool = c(0.85, 0.75, 0.65, 0.63, 0.57, 0.5, 0.43, 0.2, 0)

security <- -diff(c(1, Security))
milk <- -diff(c(1, Milk))
wool <- -diff(c(1, Wool))


#Set up for simulation
#args: reps, n, labels, models, probs1, probs2

reps <- 1000
n <- 50
models <- c("wool", "security")
probs1 <- wool
probs2 <- security

ks.sim <- function(n, models,  probs1, probs2, reps=1000) {

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
  
  #Return results
  list(x, D.max, cum.totals)
  
}

x <- ks.sim(30, reps=100000, models=c("wool", "security"), probs1=wool, probs2=security)
x[[1]]
