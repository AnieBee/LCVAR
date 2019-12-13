calculateB <- function(Covariates, K, nDepVar, A, Sigma, N, PersPDiffStart, PersEnd, X, Y, Lags, FZY, qqq, B)
{
    if(Covariates == "equal-within-clusters"){# B(0): equal-within-clusters
        for(j in 1:K)
        {
            Bnum = 0
            Bdenom = 0
            Delta = cbind(diag(nDepVar), (-1) * A[ , 1:(nDepVar * Lags[j]), j])
            NewCovaMat = t(Delta) %*% ginv(Sigma[, , j]) %*% Delta
            for(i in 1:N)
            {
                Bn = 0
                Bd = 0
                for(trunner in (PersPDiffStart[[Lags[j]]][i]):(PersEnd[i]) )
                { 
                    XtildaKron <- kronecker(t(X[ , (trunner):(trunner - Lags[j]), drop = FALSE]),
                                            diag(nDepVar)) # drop = FALSE is needed here in case of q = 1
                    Bd = Bd + ( t(XtildaKron) %*% NewCovaMat %*% XtildaKron ) 
                    Bn = Bn + ( t(XtildaKron) %*% NewCovaMat %*% 
                                     as.vector(Y[ , (trunner):(trunner - Lags[j]), drop = FALSE]) )
                }
                Bnum = Bnum + (FZY[ i, j] * Bn)
                Bdenom = Bdenom + (FZY[ i, j] * Bd)
            }
            BasVec = ginv(Bdenom) %*% Bnum
            B[, , j] = matrix(BasVec, nrow = nDepVar, ncol = qqq, byrow = FALSE)
        }
    }# else
    # {# Constraints:
    #     if(Covariates == "equal-across-clusters")
    #     { # B(0): equal-accros-clusters 
    #         ####The  code below is unchecked#####
    #        
    #          Bnum = 0
    #         Bdenom = 0
    #         for(j in 1:K)
    #         {
    #             Delta = cbind(diag(nDepVar), (-1) * A[ , 1:(nDepVar * Lags[j]), j])
    #             NewCovaMat = t(Delta) %*% ginv(Sigma[, , j]) %*% Delta
    #             for(i in 1:N)
    #             {
    #                 Bn = 0
    #                 Bd = 0
    #                 for(trunner in (PersPDiffStart[[Lags[j]]][i]):(PersEnd[i]) )
    #                 { #check this is correct in runner 
    #                     XtildaKron = kronecker(t(X[ , (trunner):(trunner - Lags[j]), drop = FALSE]), diag(nDepVar)) # drop = FALSE is needed here in case of q = 1...possibly needed in other matrices as well to prepare for univariate case
    #                     Bd = Bd + (t(XtildaKron) %*% NewCovaMat %*% XtildaKron) 
    #                     Bn = Bn + (t(XtildaKron) %*% NewCovaMat %*% as.vector(Y[ , (trunner):(trunner - Lags[j]) ]))
    #                 }
    #                 Bnum = Bnum + (FZY[ i, j] * Bn)
    #                 Bdenom = Bdenom + (FZY[ i, j] * Bd)
    #             }
    #         }
    #         BasVec = ginv(Bdenom) %*% Bnum
    #         B[, , 1] = matrix(BasVec, nrow = nDepVar, ncol = qqq, byrow = FALSE)
    #     } # else
    #     #{ # B(0): "individual-specific"
    #     #}
    # }# End of B(0) Calculation 
    
    invisible(B)
    
}

