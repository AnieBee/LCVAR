createX <- function(YLength, xFactor, Data, xContinuous)
{
    # Create a single intercept, is needed when xFactor is not null becasue in Dum the first dummy of every single factor is removed
    # At least an Intercept is always included, even if no X variables are included
    X <- matrix(1, nrow = YLength) # intercept 
    colnames(X) <- "Intercept" 
    
    if(!is.null(xFactor))
    { # Create Dum containing dummies with first dummy removed for every xFactor variable
        Data[xFactor] <- lapply(Data[xFactor], as.factor) # Make every variable that is supposed to ba a factor into a factor
        Dum <- as.matrix(Data[, xFactor])
        colnames(Dum) <- names(Data)[xFactor] # save name so it does not lose col name if it contains only a single var
        Dum = dummy_cols(Dum, remove_first_dummy = TRUE) # Dumm contains dummies for all xFactor variables:
        # DumDum = dummy_cols(DumDum, remove_first_dummy = T) #Dummy = model.matrix(~factor(Data[, xFactor[r]]))
        Dum = Dum[ , -c(1:length(xFactor))] # delete xFactor variables after you added dummies, 
        # make sure dummies are always added Behind the existing variables, otherwise the above removes the wrong variables
        # qq <- dim(Dum)[2] # number of dummies
        X <- cbind(X, Dum) # combine Intercept with all dummy variables
    }
    
    if(!is.null(xContinuous))
    { # Make X0 containing continious variables and add them to X (contains intercept and possibly dummies if there are any)
        X0 <- as.matrix(Data [ , xContinuous])
        colnames(X0) <- names(Data)[xContinuous]
        X <- cbind(X, X0) # x tilda
    }

    invisible(t(as.matrix(X))) 
    
}