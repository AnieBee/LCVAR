
### Vary ###
numberVariables <- 4
numberContiniousVars <- 1
numberExtraVars <- 5 # + 1 col for true cluster membership & 2 for continious or categorical var & 1 for ID & 1 for timepoint
nPersons <- c(120)
# 1:numberVariables : endogenous Variables,
# numberVariables +1 : true cluster membership,
# numberVariables +2 : column containing ID,
# numberVariables +3 : column containing timepoint,
# numberVariables +4 : X (continious Var, categorical Var or filled with NAss)

### Vary ###
a = 0
nLags <- c(1) # index: i
nClusters <- c(2) # index: j
ClusterSize <- c(0) #index: k
Observations <- c(150)  #index: l
Distance <- c(2) # index: m
EuclideDistance <- c(.12)
Replications <- 1 #index: r
StatusOfB <- c(4) # 1 = intercept only,  2 = categorical var included, 3 = continious var included, 4  =  continious and categorical, 5 = group differences in intercept)
numbPhiOffDiag <- ((numberVariables ^ 2) - numberVariables) / 2
#VariancePhi_kqs <- .025
numberCovariates <- 4

SigmaU <- diag(numberVariables) + .5
#Sigma is the variance/covariance of the white noise series: set to identity

# Continious B (Bcoeffs) are equal across clusters
Bcoeffs <-  cbind( matrix((1:numberVariables) * .2, numberVariables, numberContiniousVars),
                   matrix((1:numberVariables) * .1, numberVariables, numberContiniousVars),
                   matrix((1:numberVariables) * -.15, numberVariables, numberContiniousVars),
                   matrix((1:numberVariables) * -.25, numberVariables, numberContiniousVars))
#as.matrix(Bcoeffs[ , n]) gives continious Bcoeffs for a certain cluster


BCoeffMatrixForComparisons <- vector("list", length(StatusOfB))

BCoeffMatrixForComparisons[[1]] <- NA
BCoeffMatrixForComparisons[[2]] <- cbind(rep(0, numberVariables), rep(2, numberVariables), rep(3, numberVariables), rep(NA, numberVariables))
BCoeffMatrixForComparisons[[3]] <- cbind(rep(0, numberVariables), rep(NA, numberVariables), rep(NA, numberVariables), as.matrix(Bcoeffs[ , 1]))
BCoeffMatrixForComparisons[[4]] <- cbind(rep(0, numberVariables), rep(2, numberVariables), rep(3, numberVariables), as.matrix(Bcoeffs[ , 1]))
BCoeffMatrixForComparisons[[5]] <- array(NA, dim = c(numberVariables, numberCovariates, max(nClusters)))
for(SomeRunnerWhatever in 1:(max(nClusters))){
    BCoeffMatrixForComparisons[[5]][, , SomeRunnerWhatever] <- cbind(rep(0, numberVariables),  seq(2, by = .5, length.out = numberVariables),
                                                                     seq(3, by = .5, length.out = numberVariables), as.matrix(Bcoeffs[ , SomeRunnerWhatever]))
}
# BCoeffMatrixForComparisons[[s]] for s == 5, BCoeffMatrixForComparisons[[5]][, , n] gives coefficient for cluster n
