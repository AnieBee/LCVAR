
determineLagOrder <- function (Lags, K, ClusterVARcoeffs, nDepVar)
    # Assumes Lags to be ordered from lowest to largest lag order
{ 
   
    LagSplits = unique(Lags)
     ## Do all clusters have same lag order anyway? then return 1:K
    if (length(LagSplits) ==  1)
    {
        invisible(1:K)
    }
    else
    {
        remainClusters = 1:K
        nClustersInSplits = table(Lags)
        orderOfCls = rep(NA, K)
        orderOfClsCounter = 1
        
        stopifnot(length(LagSplits) > 1) # Only call if there are any splits
        
        for (largerLagInSplitIdx in length(LagSplits):2) #  until to 2, because there are only length(LagSplits) -  1 many split
            # function loops over all larger lags in a split, not over all lags
        {
            absLagCoeffs = rep(NA, length(remainClusters))
            
            
            absLagCoeffsCounter = 0
            for (clustCount in remainClusters)
            {
                absLagCoeffsCounter = absLagCoeffsCounter + 1
                absLagCoeffs[absLagCoeffsCounter] = sum(abs(ClusterVARcoeffs[ , 
                                                                              ((LagSplits[largerLagInSplitIdx - 1] * nDepVar) + 1):
                                                                                  (LagSplits[largerLagInSplitIdx] * nDepVar)
                                                                              , clustCount]))
                
            }
            
            selectedClsIdx = order(absLagCoeffs, decreasing = TRUE)[1:(nClustersInSplits[largerLagInSplitIdx])] 
            # The Nsplit (number of clusters in current split) with highest lag coeffs in currently evaluated Lag(s)
            selectedCls = remainClusters[selectedClsIdx]
            # Which clusters belong in the current lag order?
            remainClusters =  remainClusters[-selectedClsIdx]
            # update the which clusters remain after current clusters are out of the selction
            orderOfCls[orderOfClsCounter:(orderOfClsCounter + nClustersInSplits[largerLagInSplitIdx] - 1)] = selectedCls
            orderOfClsCounter = orderOfClsCounter + nClustersInSplits[largerLagInSplitIdx]
            # Add currently selected clusters to the order vector
            
        }
        # Allocate remaining clusters to lowest lag 
        orderOfCls[orderOfClsCounter:(orderOfClsCounter + nClustersInSplits[1] - 1)] = remainClusters
        
        # check every cluster is allocated only once
        stopifnot(table(orderOfCls) == 1)
        
        invisible(rev(orderOfCls))
        # returned vector gives clusters from lowest lag order to highest lag order
    }
   
}