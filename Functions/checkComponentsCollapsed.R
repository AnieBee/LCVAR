checkComponentsCollapsed <- function(K, N, FZY, smallestClN, EMiteration, crisp = FALSE)
{

    resetCl = NULL
    ComponentColapsedOntoSinglePoint = which(table(factor(apply(FZY, 1, which.max), levels = as.character(1:K))) < smallestClN)
    # ComponentColapsedOntoSinglePoint is true if cluster contains less than smallestClN of people
    while (as.logical(length(ComponentColapsedOntoSinglePoint)))
    {  # If one cluster is empty or contains less than smallestClN: ressample everyone

        cat(c("\n Warning: A single/empty cluster occured in EM-iteration",
              EMiteration, ", memberships and Sigma reset \n"))
        cat(c("\n Warning: A single/empty cluster occured in EM-iteration",
              EMiteration, ", memberships and Sigma reset \n"), file = "EMwarnings.txt", append = TRUE) 
        resetCl = unique(c(resetCl, ComponentColapsedOntoSinglePoint))
        for (clust in ComponentColapsedOntoSinglePoint)
        {
            # FZY[order(FZY[ , clust], decreasing = TRUE)[1:smallestClN] , clust] = 1 + 1e-100 # the  highest posteriors in the empty cluster are set to 1
            FZY[sample(1:N, smallestClN, replace = FALSE), clust] = 1.01
        }
        FZY = t(scale(t(FZY), center = FALSE, scale = rowSums(FZY))) # Scale posteriors so they sum to 1 again
        
        if (crisp)
        {
            classification = apply(FZY, 1, which.max)
            for (indv in 1:N)
            {
                FZY[indv, ] = rep(0, K)
                FZY[indv, classification[indv]] = 1    
            }
        }

        ComponentColapsedOntoSinglePoint = which(table(factor(apply(FZY, 1, which.max), levels = as.character(1:K))) < smallestClN)
    }

        invisible(list(FZY = FZY, resetCl = resetCl,
                   iterationReset = as.logical(length(resetCl))))
}