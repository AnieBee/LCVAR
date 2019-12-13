calculateTau <- function(K, tau, FZY, N)
{
    for (j in 1:K)
    {
        tau[j] = sum(FZY[, j]) / N
        
    }
    
    invisible(tau)
}