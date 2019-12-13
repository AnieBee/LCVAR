rewrittenUpperTri = function (x, diag = FALSE, ColMax = 0, ColMin = 0) 
{ #is upper.tri base func rewritten for Var matrices
    x <- as.matrix(x)
    if (ColMax) 
    (row(x) < col(x) & col(x) < ColMax)
    else (row(x) < (col(x) - ColMin) & col(x) > ColMin)
}

# rewrittenUpperTri(Phi_k, ColMax = numberVariables + 1) # Upper Triangle of Lag 1 coefficents
# rewrittenUpperTri(Phi_k, ColMin = numberVariables) # Upper Triangle of Lag 2 coefficents

rewrittenLowerTri <- function (x, diag = FALSE, ColMax = 0, ColMin = 0) 
{
    x <- as.matrix(x)
    if (ColMax) 
        (row(x) > col(x) & col(x) < ColMax)
    else (row(x) > (col(x) - ColMin) & col(x) > ColMin)
}

# rewrittenLowerTri(Phi_k, ColMax = numberVariables + 1) # Lower Triangle of Lag 1 coefficents
# rewrittenLowerTri(Phi_k, ColMin = numberVariables) # Upper Triangle of Lag 2 coefficents