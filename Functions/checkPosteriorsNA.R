
checkPosteriorsNA <- function(FZY, K, EMiteration)
{
    if(any(is.na(FZY)))
    { # FZY = ifelse(is.na(FZY), 1/K, FZY) 
      # if posterior prob for any person is NA reset that posterior tau_{ij} to a "blind" prior: 1/k
        if(all(is.na(FZY))) stop("All posteriors are NA in EM-iteration", hh, "\n")
        cat("\n Warning: Some posteriors are NA in EM-iteration", EMiteration, ", memberships reset \n")
        cat("\n Warning: Some posteriors are NA in EM-iteration", EMiteration, ", memberships reset \n", file = "EMwarnings.txt", append = TRUE)
        FZY[which(is.na(FZY))] = 1 / K # Set those p(Y|z_{ik}) that are NA for some person to 1/K
        FZY = t(scale(t(FZY), center = FALSE, scale = rowSums(FZY))) # Scale posteriors so they sum to 1 again
    }

    stopifnot(sum(FZY) != 0) # sum(FZY should equal N)
    
    invisible(list(FZY = FZY, iterationReset = any(is.na(FZY))))
    
}