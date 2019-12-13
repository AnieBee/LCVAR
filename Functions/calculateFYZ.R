calculateFYZ <- function(K, N, FYZ, U, PersStartU, PersEndU, nDepVar, Sigma, Lags)
{
    for(j in 1:K){  
        for(i in 1:N)
        {
            FYZ[i, j] = sum(dmvnorm(x = t(U[ , (PersStartU[[Lags[j]]][i]):(PersEndU[[Lags[j]]][i]), j]), 
                                    mean = as.matrix(rep(0, nDepVar), ncol = nDepVar),
                                    sigma = as.matrix(Sigma[, , j], ncol = nDepVar) , log = TRUE)) 
            # Sum of all log probabilities (from normal pdf) of columns of U, 
            # columns of x represent time points per individual 
        }
    } 
    
    invisible(FYZ)
}