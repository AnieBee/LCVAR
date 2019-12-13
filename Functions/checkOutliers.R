checkOutliers <- function(K, FZY, EMiteration)
{
    Outliers = NULL
    for(j in 1:K) Outliers = unique(c(Outliers, which(!is.finite(FZY[ , j])))) # if posterior prob are infinity: 
    # because FYZandPrior is low (usually for all clusters), low FYZandPrior means low membership likelihood
    if(length(Outliers))
    { ## Check ##
        cat(c("\n In EM-iteration", EMiteration, " Outliers:", Outliers, ", memberships reset \n")) # Anja: this can go, find a better way to return this info
        cat(c("\n In EM-iteration", EMiteration, " Outliers:", Outliers, ", memberships reset \n"), file = "EMwarnings.txt", append = TRUE)
        FZY[Outliers, ] = rep(1 / K, K) # let outliers contribute equally to all clusters # Anja: in C this would need a for loop cause FZY is a matrix
    } 
    
    invisible(list(FZY = FZY, iterationReset = as.logical(length(Outliers))))
}