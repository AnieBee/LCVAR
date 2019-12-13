varsRootsFunc <- function (VARcoeffs, numberVariables, numberLags) 
{
    # A VAR(p)-process is stable, if its reverse characteristic polynomial has no roots in or on the complex circle.
    # all eigenvalues of the companion matrix A have modulus less than 1
    
    K <- numberVariables
    p <- numberLags
    A <- VARcoeffs
    companion <- matrix(0, nrow = K * p, ncol = K * p)
    companion[1:K, 1:(K * p)] <- A
    if (p > 1) {
        j <- 0
        for (i in (K + 1):(K * p)) {
            j <- j + 1
            companion[i, j] <- 1
        }
    }
    roots <- eigen(companion)$values
    roots <- Mod(roots)
    return(all(roots < 1))
    # returns TRUE if VAR(p) is stationary
}