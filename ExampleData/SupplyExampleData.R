# source("/home/anja/Desktop/PhD/LCVARPackage/ExampleData/DataGeneration.R")
load("/home/anja/Desktop/PhD/LCVARPackage/ExampleData/ScrambledData.RData")

## Look at true data
Phi_ks
BCoeffMatrixForComparisons[[5]]


Data = ScrambledData
ID = 6
yVars = 1:4
Time = 7
xContinious =  8 # NULL #
xFactor = 9
Covariates = "equal-within-clusters"
Clusters =  2
seme = 3
LowestLag = 1
HighestLag = 1
smallestClN = 3
ICType = "HQ"
Rand = 1
Rational = TRUE
Initialization = NULL
SigmaIncrease = 10
it = 25
Conv = 1e-06

Result = LCVARclust(Data  = ScrambledData, yVars = 1:4, Time = 7, ID = 6, xContinious = 8, xFactor = 9,
                                Covariates = "equal-within-clusters", # "equal-accros-clusters", "individual-specific"),  so far only equal-within-clusters is implemented
                                Clusters = 2, LowestLag = 1, HighestLag = 1, smallestClN = 3,
                                # Smallest allowed cluster, smallestClN is used in checkComponentsCollapsed
                                ICType = "HQ", seme = 3,
                                Rand = 2, Rational = T, Initialization = NULL,
                                SigmaIncrease = 10, it = 25, Conv = 1e-06)

TrueId = Data[PersStart, 5]
mem = InitMT$memb
mem[2, ] = mem[2, ] * 2
mem = mem[1, ] + mem [2, ]
