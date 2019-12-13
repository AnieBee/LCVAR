require(fastDummies)  # (v.1.1.0 ) for dummy_cols()
require(MASS) # (v.7.3-49) for ginv()
require(mvtnorm) # (v.1.0-6) for dmvnorm()

source("Functions/calculateA.R")
source("Functions/calculateB.R")
source("Functions/calculateBandWZero.R")
source("Functions/calculateCoefficientsForRandoAndRational.R")
source("Functions/calculateFYZ.R")
source("Functions/calculatePosterior.R")
source("Functions/calculateTau.R")
source("Functions/calculateRatio.R")
source("Functions/calculateNPara.R")
source("Functions/calculateIC.R")
source("Functions/calculateW.R")
source("Functions/calculateSigma.R")
source("Functions/calculateU.R")
source("Functions/calculateSigma.R")
source("Functions/calculateLagList.R")
source("Functions/callCalculateCoefficientsForRandoAndRational.R")
source("Functions/callEMFuncs.R")
source("Functions/checkComponentsCollapsed.R")
source("Functions/checkConvergence.R")
source("Functions/checkOutliers.R")
source("Functions/checkSingularitySigma.R")
source("Functions/constraintsOnB.R")
source("Functions/checkLikelihoodsNA.R")
source("Functions/checkPosteriorsNA.R")
source("Functions/createOutputList.R")
source("Functions/createX.R")
source("Functions/determineLagOrder.R")
source("Functions/reorderLags.R")
source("Functions/InitFuncs.R")
source("Functions/EMInit.R")
source("Functions/EMFunc.R")
source("Functions/LCVARclust.R")

load("ExampleData/ExampleData.RData")
head(Dataset)

Result = LCVARclust(Data = Dataset, yVars = 1:4, Time = 6, ID = 5, xContinious = 7, xFactor = 8,
                    Covariates = "equal-within-clusters",
                    Clusters = 2, LowestLag = 1, HighestLag = 2, smallestClN = 3,
                    ICType = "HQ", seme = 3, Rand = 2, Rational = TRUE, 
                    SigmaIncrease = 10, it = 25, Conv = 1e-06, Initialization = NULL)

Result$BestSolutionsPerCluster

