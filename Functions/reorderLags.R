reorderLags <- function(Lags, K, newOrder)
{
    newLags = Lags
    for (clustCounter in 1:K)
    {
        newLags[newOrder[clustCounter]] =  Lags[clustCounter]
        
    }
    
    invisible(newLags)
    
}