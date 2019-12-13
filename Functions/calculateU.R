calculateU <- function(K, WkNumbVersions, N, PersPDiffStart, PersEnd, U, Wk, A, Lags, nDepVar)
{

    for(j in 1:K)
    {
        # Determine runner for W based on the constraints on B -------------------------
        WkRunner = ifelse(WkNumbVersions == K, j, 1) 
        
        # calculate U ------
        Urunner = 0
        for(i in 1:N)
        {
            for(trunner in (PersPDiffStart[[Lags[j]]][i]):PersEnd[i])
            { # from Lags to T
                Urunner = Urunner + 1  
                # either A will always index the needed positions of itself or the runner in Wk does not have to change
                # If you would not index in A, you would not use all the Us you have at your disposal when lag number is smaller 
                U[ , Urunner, j] = Wk[ , trunner, WkRunner] - 
                                    ( A[ , 1:(nDepVar * Lags[j]), j] %*% as.vector(Wk[ , (trunner - 1):(trunner - Lags[j]), WkRunner]) ) 
            }

        }  
    }
    
    invisible(U)
    
}