# Define function to simulate two-sample KS test and X2 results
payne.simulate <- function(n, probs, breaks, reps=1000, model.names=NULL) {

    #Set up bin labels
    labels <- numeric(0)
    for(i in 1:(length(breaks) - 1)) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
    }

    #Set up model names
    if(!class(probs) == "list") {probs <- list(probs)}
    if(is.null(model.names)) {model.names <- names(probs)}
    if(is.null(model.names)) {model.names <- as.character(1:length(probs))}

    #Simulate
    rep.no <- rep(1:reps, each = n)
    totals <- data.table(rep.no=rep(1:reps, each=length(labels)), bin.no=rep(1:length(labels), reps), bin=rep(labels, reps))
    sim <- data.table(rep.no=rep(1:reps, each=n), id=rep(1:n, reps))
    for(i in 1:length(probs)) {
        sim[, bin.no := sample(1:length(labels), size=nrow(sim), replace=TRUE, prob=probs[[i]])]
        x <- sim[, length(id), by=list(rep.no, bin.no)]
        setnames(x, "V1", model.names[i])
        totals <- merge(totals, x, by=c("rep.no", "bin.no"), all=TRUE)
    }
    totals[is.na(totals)] <- 0

    #Summarise
    summary <- sim.summ(totals)

    #Return results
    list(totals, summary)

}


sim.sig.2 <- function(n, models,  probs1, probs2, reps=1000) {

    #Set up bin labels
    labels <- numeric(0)
    for(i in 1:(length(breaks) - 1)) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
    }

    #Simulate from each model
    rep.no <- rep(1:reps, each = n)
    sim <- data.table(rep.no=rep(1:reps, each=n), id=rep(1:n, reps))
    sim[, a := sample(1:length(labels), size=nrow(sim), replace=TRUE, prob=probs1)]
    sim[, b := sample(1:length(labels), size=nrow(sim), replace=TRUE, prob=probs2)]

    #Two-sample K-S test
    KS.fun <- function(sim) {
        x <- ks.test(sim$a, sim$b)
        x[[2]]
    }
    KS <- data.table(p=ddply(sim, .(rep.no), KS.fun))
    KS[, sig := 0]
    KS[p.V1 < 0.05, sig := 1]
    x <- sum(KS$sig)/nrow(KS)

    #Calculate totals per bin
    totals1 <- sim[, j=list(count.a=length(id)), by=list(rep.no, a)]
    totals2 <- sim[, j=list(count.b=length(id)), by=list(rep.no, b)]
    setnames(totals1, "a", "b")
    totals <- merge(totals1, totals2, by=c("rep.no", "b"), all=TRUE)

    #Tidy up results
    frame <- data.table(rep.no=rep(1:reps, each=length(labels)), b=rep(1:length(labels), reps), bin=rep(labels, reps))
    totals <- merge(frame, totals, by=c("rep.no", "b"), all=TRUE)
    totals[is.na(totals)] <- 0

    #Chi-squared
    X2.fun <- function(totals) {
        x <- chisq.test(totals[,c("count.a", "count.b")])
        x[[3]]
    }
    X2 <- data.table(p=daply(totals, .(rep.no), X2.fun))
    X2 <- X2[!is.na(p),]
    X2[, sig := 0]
    X2[p < 0.05, sig := 1]
    y <- sum(X2$sig)/nrow(X2)

    #Return results
    c(x,y)

}

# Define function to simulate 1-sample K-S and X2 results
sim.sig.1 <- function(n, probs, ref.models, breaks, reps=1000, model.names=NULL) {

    #Set up bin labels
    labels <- numeric(0)
    for(i in 1:(length(breaks) - 1)) {
        labels[i] <- paste(breaks[i], breaks[i + 1], sep="-") #sets bin labels
    }

    #Set up reference model names
    if(!class(ref.models) == "list") {ref.models <- list(ref.models)}
    if(is.null(model.names)) {model.names <- names(ref.models)}
    if(is.null(model.names)) {model.names <- as.character(1:length(ref.models))}
    models <- list()
    for(i in 1:length(ref.models)) {
        x <- cbind(ref.models[[i]], 1:length(ref.models[[i]]))
        models[[i]] <- ecdf(rep(x[,2], x[,1]*100000))
    }

    #Simulate
    rep.no <- rep(1:reps, each = n)
    sim <- data.table(rep.no=rep(1:reps, each=n), id=rep(1:n, reps))
    sim[, bin.no := sample(1:length(labels), size=nrow(sim), replace=TRUE, prob=probs)]

    #One-sample K-S test
    KS1.fun <- function(sim) {
        z <- ks.test(sim$bin.no, models[[i]])
        z[[2]]
    }

    x <- list()
    for(i in 1:length(ref.models)) {
        KS <- data.table(p=ddply(sim, .(rep.no), KS1.fun))
        KS[, sig := 0]
        KS[p.V1 < 0.05, sig := 1]
        x[[i]] <- sum(KS$sig)/nrow(KS)
    }
    x
}


