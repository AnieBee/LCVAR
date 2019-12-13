calculateRatio <- function(temp, lik, EMiteration, Conv, L2iterationReset)
{
    if ((temp < lik) & (EMiteration > 3) & ( ! L2iterationReset))
    { # Likelihood decreased indicating precision has been reached (not due to having been reset in this iteration)
        invisible(Conv)
        
    }else
    {
        invisible(abs((temp - lik) / lik) )
    }
}