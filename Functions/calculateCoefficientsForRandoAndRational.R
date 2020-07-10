calculateCoefficientsForRandoAndRational <- function(Covariates, K, N, nDepVar, qqq, Lag, PersEnd, PersStart, Y, X, PersPDiffStart)
{
    DimensionsBasedonConstraints = constraintsOnB(Covariates, K, N)
    BIndividual = array(NA, dim = c(nDepVar, qqq, N))
    CoefficientsForRandoAndRational = array(NA, dim = c((nDepVar * nDepVar * Lag) +
                                                            (DimensionsBasedonConstraints$ClusterOnB * (nDepVar * qqq)),
                                                        N))
    WIndividual = array(NA, dim = c(nDepVar, PersEnd[N]))
    for(i in 1:N)
    {
        BIndividual[ , , i] = (Y[ , PersStart[i]:PersEnd[i]] %*% t(X[ , PersStart[i]:PersEnd[i], drop = FALSE])) %*%
            ginv(X[ , PersStart[i]:PersEnd[i], drop = FALSE] %*% t(X[ , PersStart[i]:PersEnd[i], drop = FALSE]))
        WIndividual[ , (PersStart[i]):(PersEnd[i])] = Y[ , (PersStart[i]):(PersEnd[i])] - (BIndividual[ , , i] %*% X[ , (PersStart[i]):(PersEnd[i]), drop = FALSE])
    }
    for(i in 1:N)
    {
        AKn = 0
        AKd = 0
        for(trunner in PersPDiffStart[i]:PersEnd[i])
        {
            AKn = AKn + WIndividual[ , trunner] %*% t(as.vector(WIndividual[ , (trunner - 1):(trunner - Lag)])) # W[, (trunner-1):(trunner - Lag)] = Z_{it}
            AKd = AKd + as.vector(WIndividual[ , (trunner - 1):(trunner - Lag)]) %*% t(as.vector(WIndividual[ , (trunner - 1):(trunner - Lag)]))
        }
        CoefficientsForRandoAndRational[1:(nDepVar * nDepVar * Lag), i] = as.vector(AKn %*% ginv(AKd)) # AIndividual # gets read in left to right
    }
    if(DimensionsBasedonConstraints$ClusterOnB)
    {# Add B to these coefficients on which the initial clustering solutions are built
        CoefficientsForRandoAndRational[(1 + (nDepVar * nDepVar * Lag)):((nDepVar * qqq) + (nDepVar * nDepVar * Lag)),  ] = as.vector(BIndividual)
        # gets read in top to bottom
    }

    invisible(CoefficientsForRandoAndRational)

}
