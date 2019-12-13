callCalculateCoefficientsForRandoAndRational <- function(Covariates, K, N, nDepVar, qqq, HighestLag, LowestLag, PersEnd, PersStart, Y, X, PersPDiffStart)
{
    
    CoeffsForRandoAndRationalList = vector(mode = "list", HighestLag) # list elements of 1 to LowestLag are empty
    for(lagCounter in LowestLag:HighestLag)
    {
        CoeffsForRandoAndRationalList[[lagCounter]] = 
            calculateCoefficientsForRandoAndRational(Covariates = Covariates, K = K, N = N,
                                                     nDepVar = nDepVar, qqq = qqq, Lag = lagCounter, 
                                                     PersEnd = PersEnd, PersStart = PersStart, Y = Y,
                                                     X = X, PersPDiffStart = PersPDiffStart[[lagCounter]])
    }
    
    invisible(CoeffsForRandoAndRationalList)
    
}