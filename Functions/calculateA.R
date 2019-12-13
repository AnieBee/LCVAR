calculateA <- function(K, WkNumbVersions, N, Wk, PersPDiffStart, PersEnd, Lags, FZY, A, nDepVar)
{
    for (j in 1:K)
    {
        WkRunner = ifelse(WkNumbVersions == K, j, 1) 
        
        AKnum = 0 # Sum of individal sums for A weighted by memb (tau), within j 
        AKdenom = 0 
        for (i in 1:N)
        {
            AKn = 0 # individual sum for A, within j
            AKd = 0
            for (trunner in (PersPDiffStart[[Lags[j]]][i]):PersEnd[i])
            {
                AKn = AKn + Wk[, trunner, WkRunner] %*% t( as.vector(Wk[ , (trunner - 1):(trunner - Lags[j]), WkRunner]) ) # W[, (trunner-1):(trunner - Lags[j])] = Z_{it}
                AKd = AKd + as.vector(Wk[ , (trunner - 1):(trunner - Lags[j]), WkRunner]) %*% 
                    t(as.vector(Wk[ , (trunner - 1):(trunner - Lags[j]), WkRunner]))
            }
            AKnum = AKnum + (FZY[ i, j]*AKn)
            AKdenom = AKdenom + (FZY[ i, j]*AKd)
        }
        A[ , 1:(nDepVar * Lags[j]), j] = AKnum%*%ginv(AKdenom) 
    }
    
    invisible(A)
    
}