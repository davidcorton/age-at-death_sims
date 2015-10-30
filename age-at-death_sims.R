library(data.table) # speeds things up!
library(plyr) # so we can use the split-apply-combine strategy with daply/ddply
library(dgof) # upgrades ks.test to deal with one-sample tests against step functions
library(devtools) # allows loading direct from github
install_github("davidcorton/archSeries")
library(archSeries) # for plotting functions
source("age-at-death_functions.R")

# Define breaks (just for axis labels)
breaks <- c(0,2,6,12,24,36,48,72,96,120)

# Define models in terms of survivorship [NEED TO ADD MORE]
security = c(0.904, 0.904, 0.645, 0.38, 0.25, 0.239, 0.171, 0.118, 0)
milk = c(0.47, 0.42, 0.39, 0.35, 0.28, 0.23, 0.18, 0.1, 0)
wool = c(0.85, 0.75, 0.65, 0.63, 0.57, 0.5, 0.43, 0.2, 0)
bmilk = c(0.83, 0.5, 0.36, 0.18, 0.18, 0.06, 0.01, 0.01, 0)

# Convert to mortality
security <- -diff(c(1, security))
milk <- -diff(c(1, milk))
wool <- -diff(c(1, wool))
bmilk <- -diff(c(1, bmilk))

# Demonstrate concept and confidence intervals
sizes <- c(100, 75, 50, 40, 30, 20))
model.list <- list(security=security, milk=milk, wool=wool, bmilk=bmilk)
models.fun <- function(n, model.list, breaks, reps) {
    sim <- payne.simulate(n, model.list, breaks, 1000)
    x <- surv.convert(sim)
    poly.chron(x)
    legend("bottomleft", legend = paste("n =", n), bty = "n")
}
lapply(sizes, models.fun, model.list=model.list, breaks=breaks, reps=1000)

#Simulate one sample and compare against all models
sizes <- seq(2, 100, 2)

##Security
model.list <- list(security=security, milk=milk, wool=wool, bmilk=bmilk)
security.tests <- 100 * (1 - unlist(sapply(sizes, sim.sig.1, probs=security, ref.models=model.list, breaks=breaks, reps=1000)))
plot(rep(sizes, length(model.list)), security.tests, main="Security versus...", ylab = "(False) negative rate (%)", xlab = "Sample size", type="n")
col.list = c("grey80", "darkblue", "darkgreen", "darkgoldenrod")
for(i in 1:length(model.list)) {
    points(sizes, security.tests[seq(i, length(security.tests), length(model.list))], col=col.list[i], pch=19)
}
legend("right", legend=names(model.list), col=col.list[1:length(model.list)], pch=19)

##Milk
milk.tests <- 100 * (1 - unlist(sapply(sizes, sim.sig.1, probs=milk, ref.models=model.list, breaks=breaks, reps=1000)))
plot(rep(sizes, length(model.list)), milk.tests, main="Milk versus...", ylab = "(False) negative rate (%)", xlab = "Sample size", type="n")
col.list = c("darkred", "grey80", "darkgreen", "darkgoldenrod")
for(i in 1:length(model.list)) {
    points(sizes, milk.tests[seq(i, length(milk.tests), length(model.list))], col=col.list[i], pch=19)
}
legend("right", legend=names(model.list), col=col.list[1:length(model.list)], pch=19)

##Wool
wool.tests <- 100 * (1 - unlist(sapply(sizes, sim.sig.1, probs=wool, ref.models=model.list, breaks=breaks, reps=1000)))
plot(rep(sizes, length(model.list)), wool.tests, main="Wool versus...", ylab = "(False) negative rate (%)", xlab = "Sample size", type="n")
col.list = c("darkred", "darkblue", "grey80", "darkgoldenrod")
for(i in 1:length(model.list)) {
    points(sizes, wool.tests[seq(i, length(wool.tests), length(model.list))], col=col.list[i], pch=19)
}
legend("right", legend=names(model.list), col=col.list[1:length(model.list)], pch=19)

##Type B Milk
bmilk.tests <- 100 * (1 - unlist(sapply(sizes, sim.sig.1, probs=bmilk, ref.models=model.list, breaks=breaks, reps=1000)))
plot(rep(sizes, length(model.list)), bmilk.tests, main="Type B Milk versus...", ylab = "(False) negative rate (%)", xlab = "Sample size", type="n")
col.list = c("darkred", "darkblue", "darkgreen", "grey80")
for(i in 1:length(model.list)) {
    points(sizes, bmilk.tests[seq(i, length(bmilk.tests), length(model.list))], col=col.list[i], pch=19)
}
legend("right", legend=names(model.list), col=col.list[1:length(model.list)], pch=19)

#Close look at false positives
plot(sizes, 100 - security.tests[seq(1, length(security.tests), 4)], main="False positives", ylab = "False positive rate (%)", xlab = "Sample size", col="darkred", pch=19)
abline(h=5)
points(sizes, 100 - milk.tests[seq(2, length(milk.tests), 4)], col="darkblue", pch=19)
points(sizes, 100 - wool.tests[seq(3, length(wool.tests), 4)], col="darkgreen", pch=19)
points(sizes, 100 - bmilk.tests[seq(4, length(bmilk.tests), 4)], col="darkgoldenrod", pch=19)
abline(v=31.5)
### What the absolute fuck is this all about??

# Testing for false negatives by simulating two samples from two different models
sizes <- seq(10, 150, 5)

false.negs <- 100 * (1 - sapply(sizes, sim.sig.2, models=c("wool", "security"), probs1=wool, probs2=security, reps=1000))
points(sizes, false.negs[2,], ylab = "False negative rate (%)", xlab = "Sample size", col="darkblue", pch=19)
legend("topright", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.negs <- 100 * (1 - sapply(sizes, sim.sig.2, models=c("milk", "security"), probs1=milk, probs2=security, reps=1000))
plot(sizes, false.negs[1,], main="Milk vs. Security", ylab = "False negative rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.negs[2,], ylab = "False negative rate (%)", xlab = "Sample size", col="darkblue", pch=19)
legend("topright", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.negs <- 100 * (1 - sapply(sizes, sim.sig, models=c("milk", "wool"), probs1=milk, probs2=wool, reps=1000))
plot(sizes, false.negs[1,], main="Milk vs. Wool", ylab = "False negative rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.negs[2,], ylab = "False negative rate (%)", xlab = "Sample size", col="darkblue", pch=19)
legend("topright", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

# Testing for false positives
false.pos <- 100 * sapply(sizes, sim.sig, models=c("wool", "wool"), probs1=wool, probs2=wool, reps=1000)
plot(sizes, false.pos[1,], main="Wool vs. Wool", ylim=c(0,6.5), ylab = "False positive rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.pos[2,], ylab = "False positive rate (%)", xlab = "Sample size", col="darkblue", pch=19)
abline(h=5)
legend("topleft", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.pos <- 100 * sapply(sizes, sim.sig, models=c("milk", "milk"), probs1=milk, probs2=milk, reps=1000)
plot(sizes, false.pos[1,], main="Milk vs. Milk", ylim=c(0,6.5), ylab = "False positive rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.pos[2,], ylab = "False positive rate (%)", xlab = "Sample size", col="darkblue", pch=19)
abline(h=5)
legend("topleft", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")

false.pos <- 100 * sapply(sizes, sim.sig, models=c("security", "security"), probs1=security, probs2=security, reps=10000)
plot(sizes, false.pos[1,], main="Security vs. Security", ylim=c(0,6.5), ylab = "False positive rate (%)", xlab = "Sample size", col="darkred", pch=19)
points(sizes, false.pos[2,], ylab = "False positive rate (%)", xlab = "Sample size", col="darkblue", pch=19)
abline(h=5)
legend("topleft", pch=19, col=c("darkred", "darkblue"), legend=c("K-S", "X2"), bty="n")


###orphans
#Convert to cumulative
cum.totals <- totals[, j=list(a=cumsum(count.a), b=cumsum(count.b)), by="rep.no"]
cum.totals[, D := (a - b) / n]
setnames(cum.totals, c("a", "b"), models)

#Kolmogorov-Smirnov
D.max <- cum.totals[, j=list(Dmax=max(abs(D))), by="rep.no"]
cut.off <- 1.36 * sqrt((n*2)/n^2)
D.max[, sig := 0]
D.max[Dmax > cut.off, sig := 1]
x <- sum(D.max$sig)/reps

