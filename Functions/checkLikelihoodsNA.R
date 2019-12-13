checkLikelihoodsNA <- function(FYZ, EMiteration)
{
    iterationReset = FALSE
    stopifnot(sum(FYZ, na.rm = TRUE) != 0) # check for underflow, if all posteriors are zero
    if(any(is.na(FYZ)))
    { # Check not a single posterior is NA
        if(all(is.na(FYZ)))
        {
            cat("All likelihoods are NA in EM-iteration", EMiteration, "\n", file = "EMwarnings.txt", append = TRUE)
            stop("All likelihoods are NA in EM-iteration", EMiteration, "\n")
        } 
        cat("\n Warning: Some likelihoods NA in EM-iteration", EMiteration, ", likelihoods reset \n")
        cat("\n Warning: Some likelihoods NA in EM-iteration", EMiteration, ", likelihoods reset \n", file = "EMwarnings.txt", append = TRUE)
        NAIndexFYZ = which(is.na(FYZ))
        FYZ[NAIndexFYZ] = mean(FYZ, na.rm=TRUE) # Set those p(Y|z_{ik}) that are NA for some person to the overall mean of all P(Y|Z_{ik})
        iterationReset = TRUE
    }
    # medie <- colMeans(FYZ)
    # for (j in 1:K) if (is.na(medie[j])) FYZ[ ,j] <- 1e-100 # if FYZ has remaining missings, set P(Y|Z_{ik}) to a low value for every person for the cluster in question. Can only happen if
    # # mean(FYZ, na.rm=TRUE) is giving na (so all are na)
    
    invisible(list(FYZ = FYZ, iterationReset = iterationReset))
}