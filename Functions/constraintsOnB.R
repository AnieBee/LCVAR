constraintsOnB <- function(Covariates, K, N){ 
    
    if(Covariates == "equal-within-clusters")
    { # Same within clusters, dim(B)[3] = K
        BNumbVersions = K
        WkNumbVersions = K
        ClusterOnB = 1
    }else
    {# if B is constraint, Wk is same for all clusters and B is same for all clusters  dim(B)[3] = 1
        WkNumbVersions = 1
        ClusterOnB = 0
        if(Covariates == "equal-across-clusters")
        { 
            BNumbVersions = 1
        }else
        { # it ends here if "individual-specific", different for every individual dim(B)[3] = N

        }
    }
    
    invisible(list(BNumbVersions = BNumbVersions, WkNumbVersions = WkNumbVersions,
                   ClusterOnB = ClusterOnB))
    
}