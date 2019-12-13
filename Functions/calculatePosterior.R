calculatePosterior <- function(N, K, tau, FYZ, lowest.Likelihood, EMiteration, FZY, UseFZY)
{
    FYZandPrior = matrix(0, N, K)
    PosteriorDenom = matrix(0, N)
    iterationReset = FALSE  
    
    for(clust in 1:K) FYZandPrior[, clust] = log(tau[clust]) + FYZ[, clust]
    for(indv in 1:N)
    {
        AexpTrick = max(FYZandPrior[indv, ]) # A in the exp trick
        PosteriorDenom[indv, ] = AexpTrick + log(sum(exp(FYZandPrior[indv, ] - AexpTrick)))
    }
    
    if (UseFZY == FALSE)
    { # If liklelihood is used after function call, PosteriorDenom is replaced with high value,
        # if FZY are used they are checked later for infinite values so you can let infinite values be 
        if(any(is.infinite(PosteriorDenom)))
        { # check: no -Inf likelhioods, if there are replace with lowest.Likelihood specified at begining of EMfunc
            cat("\n Infinite likelihood in EM-iteration", EMiteration, ", likelihood reset \n")  
            cat("\n Infinite likelihood in EM-iteration", EMiteration, ", likelihood reset \n", file = "EMwarnings.txt", append = TRUE) 
            PosteriorDenom[which(is.infinite(PosteriorDenom))] = lowest.Likelihood 
            # PosteriorDenom = ifelse(is.finite(PosteriorDenom), PosteriorDenom, lowest.Likelihood) 
            # check: no -Inf likelhioods, if there are replace with lowest.Likelihood specified at begining of EMfunc
            iterationReset = TRUE  
        }
    }
    
    for(clust in 1:K)
    { # calculate FZY = posterior(pi_{ik})
        FZY[ , clust] = exp(FYZandPrior[ , clust] - PosteriorDenom) # calculate log metric posterior and take exp() to have real posterior in FZY
    }
    
    invisible(list(FZY = FZY, logLikelihood = sum(PosteriorDenom), iterationReset = iterationReset))
}

