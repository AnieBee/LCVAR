checkSingularitySigma <- function(nDepVar, K, Sigma, EMiteration)
{
    iterationReset = FALSE
    if(nDepVar > 1)
    {
        for(j in 1:K)
        {
            if(det(Sigma[, , j]) < 1.0e-200)
            { # Make sure S is invertible (i.e. component is not collapsed onto a single data point)
                diag(Sigma[, , j]) = diag(Sigma[, , j]) + 0.01
                iterationReset = TRUE
                cat("\n Warning: Sigma reset in EM-iteration", EMiteration, ", \n")
                cat("\n Warning: Sigma reset in EM-iteration", EMiteration, ", \n", file = "EMwarnings.txt", append = TRUE)
                
            } 
        }
    }
    
    invisible(list(Sigma = Sigma, iterationReset = iterationReset))
}