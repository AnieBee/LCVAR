selectA <- function(nDepVar, Lags, K, AforP, ICLagOrder, keepTrackLagOrder)
    # Create A matrix that contains best-fitting A of every cluster, indicated by ICLagOrder
    # Update keepTrackLagOrder according to changes in A 
{
    VARMatrix = array(0, dim = c(nDepVar, nDepVar * Lags, K))
    selectedLagOrder = array(0, dim = c(K))
    
    for (j in 1:K)
    {
        selectedLagOrder[j] = which.min(ICLagOrder[j, ])
        VARMatrix[ , , j] = AforP[ , , j, selectedLagOrder[j]] # ANJA: is min right here or do you need max
    }

    # compute selectedLagOrder to KeepTrackLagOrder$currentLagOrder and after update 
    # KeepTrackLagOrder$currentLagOrder to be selectedLagOrder
    keepTrackLagOrder$downwardChanges = which(keepTrackLagOrder$currentLagOrder > selectedLagOrder)
    keepTrackLagOrder$upwardChanges = which(keepTrackLagOrder$currentLagOrder < selectedLagOrder)
    keepTrackLagOrder$currentLagOrder = selectedLagOrder
    
    invisible(list(A = VARMatrix, keepTrackLagOrder = keepTrackLagOrder))
    
}