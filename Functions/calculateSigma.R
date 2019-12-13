calculateSigma <- function(K, N, FZY, U, PersStartU, PersEndU, Tni, Sigma, Lags)
{
    for(j in 1:K)
    {
        Snum = 0
        Sdenom = 0
        for(i in 1:N)
        { # there is a runner inside U for every individual i so the U multiplication result can be weighted by the posterior pi (memb)
            Snum = Snum + (FZY[ i, j] * (U[ , (PersStartU[[Lags[j]]][i]):(PersEndU[[Lags[j]]][i]), j] %*%
                                               t(U[ , (PersStartU[[Lags[j]]][i]):(PersEndU[[Lags[j]]][i]), j]) ) )
            Sdenom = Sdenom + (FZY[ i, j] * (Tni[[Lags[j]]][i]) )
        }
        Sigma[, , j] = Snum / Sdenom # Sdenom is integer
    }
    
    invisible(Sigma)
    
}