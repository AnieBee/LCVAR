calculateLagList <- function(K, HighestLag, LowestLag)
{
    
    lagList = vector("list", length = K) 
    for (lagRunner in 1:K)
    {
        lagList[[lagRunner]] = c(HighestLag:LowestLag)
    }
    lagList = expand.grid(lagList) 
    lagList = t(apply(lagList, 1, sort))
    lagList = lagList[!duplicated(lagList), , drop = FALSE]
    
    invisible(lagList[order(rowSums(lagList), decreasing = TRUE) , , drop = FALSE])
    # Order Lags from largest number of Paras to smallest
    
}






