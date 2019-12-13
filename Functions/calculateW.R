calculateW <- function(Covariates, K, Wk, Y, B, X)
{
    ### Calculate Wk based on Constraints on B --------
    if(Covariates == "equal-within-clusters"){ # Same within clusters, dim(B)[3] = K 
        for(j in 1:K)
        {
            Wk[ , , j] <- Y - B[ , , j] %*% X
        }
    }else{# if B is constraint, Wk is same for all clusters and B is same for all clusters dim(B)[3] = 1
        if(Covariates == "equal-across-clusters")
        { 
            Wk[ , , 1] <- Y - B[ , , 1] %*% X # if equal across clusters
        } # else{ # it ends here if "individual-specific", different for every individual dim(B)[3] = N 
        #   for(i in 1:N){
        #     Wk[ , (PersStart[i]):(PersEnd[i]), 1] <- Y - B[ , , i]%*%X[ , (PersStart[i]):(PersEnd[i])] # 
        #   }
        # }
    }

    invisible(Wk)
    
}