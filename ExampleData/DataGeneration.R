
source("/home/anja/Desktop/PhD/LCVARPackage/ExampleData/SupportFunctions/varsRootsFunc.R")
source("/home/anja/Desktop/PhD/LCVARPackage/ExampleData/SupportFunctions/globalVariables.R")

# library(devtools)
# install_version("portes", version = "2.1-3",  repos = "https://cran.rstudio.com") # varima.sim()
# save.image()

library(portes)
# library(MASS) #mvrnorm()
set.seed(705)


# SunnySideUp contains data for all individuals of a cluster after one another.
# To have a non-ordered dataset where individuals of the same cluster are
# not presented after one another, persons are scrambled unorderly. The result is
# ScrambledData which is the final data set for the condition that is saved in DataSets.
SunnySideUpData <- array(NA, dim = c(Observations,
                                     numberVariables + numberExtraVars, nPersons))
ScrambledData <- array(NA, dim = c(Observations * nPersons,
                                   numberVariables + numberExtraVars))
IDVariable <- c(sample(1:nPersons, nPersons, replace = FALSE))
PeopleAlreadyInThisDataSet <- 0
Phi_ks <- array(0, dim = c(numberVariables, numberVariables*max(nLags), 
                           length(nLags), # index:i
                           max(nClusters)))

### First Phi_k is determined for every cluster ###
# Lag 1 coeffs


Phi_k <- matrix(c(.4, -.5, .3, .4, -.1, .5, -.3, -.4, .1, .1 , .2, .2, .25, .15, .14, .5)
                , ncol = numberVariables*nLags, nrow = numberVariables)
Phi_ks[ , 1:(numberVariables * nLags), , 1] = Phi_k      # save Phi_kn
#Stationairy?:
stopifnot(varsRootsFunc(VARcoeffs = Phi_k[  , 1:(numberVariables * nLags)],
                                 numberVariables = numberVariables, numberLags = nLags))


Phi_kn = Phi_k + matrix(c(-.7, .7, .35, -.3, .4, -.3, .3, .3, .4, .15 , -.2, .2, -.45, .5, .4, -.5)
                        , ncol = numberVariables*nLags, nrow = numberVariables)
#Stationairy?:
Phi_ks[ , 1:(numberVariables * nLags), , 2] = Phi_kn       # save Phi_kn
stopifnot(varsRootsFunc(VARcoeffs = Phi_kn[  , 1:(numberVariables * nLags)],
                  numberVariables = numberVariables,
                  numberLags = nLags))


for(n in 1:(nClusters))
{
    nPerClust <- matrix((1 - ClusterSize) / (nClusters - as.numeric(ClusterSize != 0)), nrow = nClusters, ncol = 1)
    # nPerClust is vector of length nClusters[j]
    # contains the remaining percentages distributed equally
    if(ClusterSize)
    { # if first cluster has different percentage, replace equal precentage with it
        nPerClust[1] = ClusterSize
    }
    nPerClust = (nPerClust * nPersons)
    
    for(q in 1:nPerClust[n]){  # n gives the number of current evaluated cluster
       
        Phi_kq <- array(Phi_ks[ , 1:(numberVariables * nLags), , n],
                        dim = c(numberVariables, numberVariables, nLags))
        
        
        Data <- varima.sim(phi = Phi_kq, theta = NULL,  n = (Observations),
                           sigma = SigmaU) # Zero mean time series
        # Data is array[rows hold numberVariables, col hold Observations]
        
      
        ## what is the mean of the series?
        apply(Data, 2, mean)
        Data = apply(Data, 2, function(colInput) colInput - mean(colInput))   # Center at zero cause varima.sim is a nightmare when it comes to means of a series
        apply(Data, 2, mean)
        Data <- cbind(Data, rep(n, (Observations))) # Add column containing true cluster membership
        Data <- cbind(Data, rep(IDVariable[PeopleAlreadyInThisDataSet + q],  (Observations))) # Add column containing ID
        Data <- cbind(Data, 1:(Observations)) # Add column containing timepoint
        
{
            # Continious Var and Categorical Var, differences between clusters
            X1 <- rnorm( (Observations), mean = 20, sd = 20)
            Data[ , 1:numberVariables] <-  Data[, 1:numberVariables] + t(as.matrix(Bcoeffs[ , n]) %*% X1)
            X2 <- rep(c(1:3), (Observations) / 3) # make time variable: morning, midday, night: 1, 2, 3
            if((Observations) %% 3){
                X2 <- c(X2, c(1:3)[1:((Observations) %% 3)]) # fill with modulus many
            }
            Data[which(X2 == 2), 1] <-  Data[ which(X2 == 2), 1] + 2 # rnorm(1, mean = .7, sd = .05)
            Data[which(X2 == 2), 2] <-  Data[ which(X2 == 2), 2] + 2.5 # rnorm(1, mean = .7, sd = .05)
            Data[which(X2 == 2), 3] <-  Data[ which(X2 == 2), 3] + 3 # rnorm(1, mean = .7, sd = .05)
            Data[which(X2 == 2), 4] <-  Data[ which(X2 == 2), 4] + 3.5 # rnorm(1, mean = .7, sd = .05)
            
            Data[which(X2 == 3), 1] <-  Data[ which(X2 == 3), 1] + 3 # rnorm(1, mean = .7, sd = .05)
            Data[which(X2 == 3), 2] <-  Data[ which(X2 == 3), 2] + 3.5 # rnorm(1, mean = .7, sd = .05)
            Data[which(X2 == 3), 3] <-  Data[ which(X2 == 3), 3] + 4 # rnorm(1, mean = .7, sd = .05)
            Data[which(X2 == 3), 4] <-  Data[ which(X2 == 3), 4] + 4.5 # rnorm(1, mean = .7, sd = .05)
        }
        
        Data <-  cbind(Data, X1, X2)
        SunnySideUpData[ , , PeopleAlreadyInThisDataSet + q] <- Data
    }
    PeopleAlreadyInThisDataSet <- PeopleAlreadyInThisDataSet + q
}

# Scramble Sunny Side up and save in DataSets so people appear unordered in the final dataset
# Scrambled Data is ordered by ID but clusters are scrabled up, do not appear next to one another
ScrambledData <- apply(SunnySideUpData[ , , c(order(SunnySideUpData[1, numberVariables + 2, ])) ], 2, cbind)


## Look at true data
Phi_ks
BCoeffMatrixForComparisons[[5]]

save(ScrambledData, Phi_ks, BCoeffMatrixForComparisons,
     file = "/home/anja/Desktop/PhD/LCVARPackage/ExampleData/ScrambledData.RData")

