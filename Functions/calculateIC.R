calculateIC <- function(ICType, Sigma, Lags, nDepVar, K, N, FZY, Tni, tau)
    # clusterTimePoints = T - k but extended to the cluster case
    # Calculate IC for all clusters for a certain lag number (YWResidualCovariance, clusterTimePoints and LagsForP)
{

    ParasCl = as.vector(rep(0, K))
    clTimepoints = as.vector(rep(0, K))
    penaltyTerm = as.vector(rep(0, K))
    clIC = rep(0, K)
    
    for(j in 1:K)
    {
        # ParasCl[j] = calculateNPara(Lags = Lags[j], nDepVar = nDepVar, K = 1,
        #                             BnumbVersions = ifelse(BnumbVersions == 1,
        #                                                           tau[j],
        #                                                           1),
        #                             ncovariates = ncovariates)  + ((K - 1) / K) # includes all paras except tau, so + ((K - 1) / K)
        ParasCl[j] = Lags[j] * (nDepVar * nDepVar)
        for (i in 1:N)
        {
            clTimepoints[j] = clTimepoints[j] + (FZY[ i, j] * (Tni[[Lags[j]]][i]))
        }
        
        penaltyTerm[j] = switch(ICType, 
                                 "HQ" = log(log(clTimepoints[j])), # is called HQ(n) in VARselect from vars package
                                 "SC" = log(clTimepoints[j]),
                                 "AIC" = 1)
        
        ### Use ParasCl, clTimepoints and penaltyTerm to calculate clIC ###
        clIC[j] = tau[j] * ( log(det(Sigma[ , , j])) + ((2 * ParasCl[j] * penaltyTerm[j]) / clTimepoints[j]) )
        
    }
    ### Check: if any clIC is negative ###
    # clIC[which(clIC < 0)] = NA 
    if (any(clIC < 0))
    {
        cat(c("\n", "Negative clI:", clIC, "\n"))
        cat(c("\n", "Negative clI:", clIC, "\n"), file = "EMwarnings.txt", append = TRUE) 
    }
    
    invisible(sum(clIC))
    
}