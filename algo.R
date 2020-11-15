## BH procedure 
summary_BH <- function(pvals, H,
                       alphas = seq(0.01, 0.3, 0.01)){
    n <- length(pvals)
    results <- sapply(alphas, function(alpha){
        khat <- max(c(0, which(sort(pvals) <= alpha * (1:n) / n)))
        alpha <- alpha * khat / n
        nfrej <- sum(pvals[!H] < alpha, na.rm = TRUE)
        ntrej <- sum(pvals[H] < alpha, na.rm = TRUE)
        return(c(nfrej, ntrej))
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

## Storey's BH Procedure 
summary_storey <- function(pvals, H, thr = 0.5,
                           alphas = seq(0.01, 0.3, 0.01)){
    pi0 <- min(1, mean(pvals > thr) / (1 - thr))
    pvals[pvals > thr] = Inf
    n <- length(pvals)
    results <- sapply(alphas, function(alpha){
        khat <- max(c(0, which(sort(pvals) <= alpha / pi0 * (1:n) / n)))
        alpha <- alpha * khat / n
        nfrej <- sum(pvals[!H] < alpha, na.rm = TRUE)
        ntrej <- sum(pvals[H] < alpha, na.rm = TRUE)
        return(c(nfrej, ntrej))        
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

## AdaPT 
summary_adapt <- function(adapt, H, pvals){
    results <- apply(adapt$s, 2, function(s){
        tmp <- (pvals <= s)
        nfrej <- sum(tmp[!H], na.rm = TRUE)
        ntrej <- sum(tmp[H], na.rm = TRUE)
        return(c(nfrej, ntrej))
    })
    nfrej <- as.numeric(results[1, ])
    ntrej <- as.numeric(results[2, ])
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H),1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}