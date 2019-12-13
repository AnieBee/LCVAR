calculateNPara <- function(Lags, nDepVar, K, BnumbVersions, ncovariates)
{
    parasA = sum(Lags) * (nDepVar * nDepVar) # number of parameters estimated in A
    parasSigma = K * (nDepVar * (nDepVar + 1) / 2)# number of parameters estimaed in Sigma
    parasB = BnumbVersions * ncovariates * nDepVar
    nPara = (K - 1) + parasB + parasA + parasSigma # Total number of Paramters for taus, Bs, As and Sigmas 
    # Number of Parameters needed for B depend on (BNumbVersions)
    
    invisible(nPara)
    
}