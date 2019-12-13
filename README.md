LCVARclust Example
================
Anja Franziska Ernst </br> a.f.ernst\[at\]rug.nl </br>
December 12, 2019 </br>

</br> </br> This document illustrates how to fit a latent class vector-autoregressive model to a time series using the LCVARclust R function.

</br> For technical details see:

Ernst, A. F., Albers, C. J., Jeronimus, B. F. & Timmerman, M. E. (2020). Inter-individual differences in multivariate time series: Latent class vector-autoregressive modelling. *European Journal of Psychological Assessment*.

</br> Please report bugs to:
a.f.ernst\[at\]rug.nl.

</br> </br> </br>

Load the required packages
--------------------------

``` r
require(fastDummies) # (v.1.1.0) for dummy_cols()
require(MASS) # (v.7.3-49) for ginv()
require(mvtnorm) # (v.1.0-6) for dmvnorm()
```

Source all functions that are required to run LCVARclust
--------------------------------------------------------

``` r
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
```

Load the data
-------------

As an example we will analyze a generated data set "Dataset". Make sure no missing values are included in your dataframe, the LCVARclust function does not allow any missing values.

``` r
load("ExampleData/ExampleData.RData")
head(Dataset)
```

    ##         Item1    Item2     Item3     Item4 Person Timepoint
    ## [1,] 2.869908 3.755102  6.026770  6.868950      1         1
    ## [2,] 3.143062 5.199664  6.533799  9.052157      1         2
    ## [3,] 1.289307 3.893883  1.542538  4.086377      1         3
    ## [4,] 2.782861 6.461169  9.160752 13.078920      1         4
    ## [5,] 4.091992 5.336253  7.747059 10.067667      1         5
    ## [6,] 5.913368 5.124228 10.431553  9.968930      1         6
    ##      ContiniousVariable CategoricalVariableTimeOfDay
    ## [1,]          16.049110                            1
    ## [2,]           8.086699                            2
    ## [3,]          -1.524683                            3
    ## [4,]          40.518504                            1
    ## [5,]          24.089006                            2
    ## [6,]          22.051702                            3

Run the LCVARclust function
---------------------------

The following arguments can be specified:

-   **Data**: The dataframe to be used.

-   **yVars**: An integer vector specifying the position of the column(s) in dataframe **Data** that contain the endogenous variables (= the VAR time series).

-   **Time**: An integer specifying the position of the column in dataframe **Data** that contains the time point.

-   **ID**: An integer specifying the position of the column in dataframe **Data** that contains the ID variable for every participant.

-   **Covariates**: Constraints on the parameters of the exogenous variable(s). So far only "equal-within-clusters" can be specified.

-   **Clusters**: An integer vector specifying the numbers of mixture components (clusters) for which LCVAR models are to be calculated.

-   **LowestLag**: An integer specifying the lowest number of VAR(*p*) lags to consider in the calculation of LCVAR models.

-   **HighestLag**: An integer specifying the higest number of VAR(*p*) lags to consider in the calculation of LCVAR models.

-   **smallestClN**: An integer specifying the lowest number of individuals allowed in a cluster. When during estimation the crisp cluster membership of a cluster indicates less than **smallestClN** individuals, the covariance matrix and the posterior probabilities of cluster membership are reset.

-   **ICType**: The information criterion used to select the ideal model for a given number of clusters across all EM-starts and lag combinations. One of c("HQ", "SC", "AIC").

-   **seme**: An integer specifying the value supplied to `set.seed()`. Using the same seed guarantees reproducibility of solutions.

-   **Rand**: The number of pseudo-random EM-starts used in fitting each possible model.

-   **Rational**: Should a rational EM-start be used as well? Accepts TRUE or FALSE.

-   **SigmaIncrease**: A numerical value specifying the value by which every element of Sigma will be increased when posterior probabilities of cluster memberships are reset.

-   **it**: An integer specifying the maximum number of EM-iterations allowed after every EM-start. After **it** EM-iterations an EM-start is forced to terminate.

-   **Conv**: A numerical value specifying the convergence criterion of the log likelihood to determine convergence of an EM-start. For details see Ernst et al. (2020) Inter-individual differences in multivariate time series: Latent class vector-autoregressive modelling.

Optional arguments:

-   **xContinuous**: An integer vector specifying the position of the column(s) in dataframe **Data** that contain the continuous exogenous variable(s).

-   **xFactor**: An integer vector specifying the position of the column(s) in dataframe **Data** that contain the categorical exogenous variable(s).

-   **Initialization**: An integer specifying the position of a column in dataframe **Data** that contains a guess at participants' cluster membership for a fixed number of clusters.

For every fixed number of clusters as specified in **Clusters**, each combination of lag orders between **LowestLag** and **HighestLag** corresponds to a different statistical model. The models associated with all possible combinations of lag orders are estimated. For each model, the algorithm uses several EM-starts based on: pseudo-random initializations (**Rand**), a k-means based rational initialization (**Rational**), a guess at cluster membership (**Initialization**), and the use of a previous solution (always used). Every EM-start leads to one solution after several EM-iterations. Solutions are either reached because the likelihood converged (**Conv**) or because the maximum number of EM-iterations has been reached (**it**). Thus several solutions are reached for every possible statistical model. The ideal statistical model for a given number of clusters across all EM-starts and all lag combinations is determined with the information criterion specified in **ICType**. As a result, for every number of clusters specified in **Clusters** there will be one solution displayed in `Result$BestSolutionsPerCluster`.

``` r
Result = LCVARclust(Data = Dataset, yVars = 1:4, Time = 6, ID = 5, 
                    Covariates = "equal-within-clusters",
                    Clusters = 2, LowestLag = 1, HighestLag = 2, smallestClN = 3,
                    ICType = "HQ", seme = 3, Rand = 2, Rational = TRUE, 
                    SigmaIncrease = 10, it = 25, Conv = 1e-06, xContinuous = 7, xFactor = 8)
```

    ## 
    ##  2 Clusters: 
    ##   Lags: 2 2 
    ## * 1   * 2   * Rational   * Previous   
    ##   Lags: 1 2 
    ## * 1   
    ##  Warning: A single/empty cluster occured in EM-iteration 1 , memberships and Sigma reset 
    ## 
    ##  Warning: A single/empty cluster occured in EM-iteration 2 , memberships and Sigma reset 
    ## 
    ##  Warning: A single/empty cluster occured in EM-iteration 3 , memberships and Sigma reset 
    ## 
    ##  Warning: A single/empty cluster occured in EM-iteration 4 , memberships and Sigma reset 
    ## 
    ##  Warning: A single/empty cluster occured in EM-iteration 5 , memberships and Sigma reset 
    ## 
    ##  EM did not converge: Iteration terminated after reset stagnation 
    ## * 2   * Rational   * Previous   
    ##   Lags: 1 1 
    ## * 1   * 2   * Rational   * Previous

Interpret the output
--------------------

`Result` contains two lists: `BestSolutionsPerCluster` which contains the ideal solution for each number of clusters, and `AllSolutions` which contains the solutions for all starts of all estimated models. The two lists are structured as follows:

-   `Result$BestSolutionsPerCluster[[a]]`: contains the ideal solution for the ath number of clusters within range **Clusters** across all lag combinations and EM-starts.
-   `Result$AllSolutions[[a]][[b]][[c]]:` contains the solution for the ath number of clusters within range **Clusters** for the bth combination of lag orders on the cth EM-start.

The output for every solution contains:

-   **Converged**: Whether this EM-start converged.

-   **A**: The VAR(p)-coefficients for all clusters. The rows give the variables, the columns the lag coefficients. Lag coefficients quantify the influence the endogenous variables (columns) have on the endogenous variables (rows) at future time points.

-   **B**: The exogenous coefficients for all clusters.

-   **EMRepetitions**: The Number of EM-iterations before the algorithm terminated.

-   **last.loglik**: The log likelihood of the this model.

-   **nPara**: The number of parameters of this model.

-   **Sigma**: The error covariance matrix.

-   **LogLikelihood**: A numeric vector showing the value of the log likelihood after every EM-iteration.

-   **EMiterationReset**: A vector of logical values indicating whether any values were reset during an EM-iteration.

-   **PosteriorProbs**: A numeric vector indicating the posterior probabilities of cluster membership of every person. The output is ordered by factor(**ID**).

-   **Lags**: The number of VAR(*p*) lags for every cluster in this solution. In our example, the best solution for two number of clusters is based on 1 lag for both clusters, we can conclude that the information criterion specified in **ICType** selected a VAR(1) model for both clusters.

-   **Classification**: The crisp cluster membership for every person in this solution. The output is ordered by factor(**ID**).

-   **IC**: The value of the information criterion that was specified in **ICType**.

-   **SC**: The value of the SC information criterion.

-   **Proportions**: The estimated mixing proportions of the clusters.

``` r
Result$BestSolutionsPerCluster
```

    ## [[1]]
    ## [[1]]$Converged
    ## [1] TRUE
    ## 
    ## [[1]]$A
    ## , , 1
    ## 
    ##              Item1      Item2       Item3       Item4
    ## Item1  0.009113468  0.2712338 0.004804005 0.002756114
    ## Item2 -0.413364606  0.3943745 0.262143989 0.298928692
    ## Item3  0.600145815 -0.3102121 0.004383167 0.279625558
    ## Item4  0.311616188 -0.4046344 0.394373760 0.551653570
    ## 
    ## , , 2
    ## 
    ##             Item1      Item2     Item3     Item4
    ## Item1  0.08479586 -0.1033597 0.1068198 0.2466591
    ## Item2 -0.21764302  0.1776732 0.1040631 0.1632604
    ## Item3  0.29650512 -0.3046279 0.1888623 0.1450859
    ## Item4  0.40965841 -0.4093675 0.1878278 0.4930011
    ## 
    ## 
    ## [[1]]$B
    ## , , 1
    ## 
    ##          Intercept CategoricalVariableTimeOfDay_2
    ## Item1  0.017391963                       1.981168
    ## Item2  0.019302515                       2.456360
    ## Item3 -0.037663140                       3.014690
    ## Item4 -0.003983069                       3.497786
    ##       CategoricalVariableTimeOfDay_3 ContiniousVariable
    ## Item1                       2.967675         0.09994573
    ## Item2                       3.459679         0.20031637
    ## Item3                       4.032440         0.30097336
    ## Item4                       4.509843         0.39992009
    ## 
    ## , , 2
    ## 
    ##           Intercept CategoricalVariableTimeOfDay_2
    ## Item1 -0.0007224385                       2.001126
    ## Item2  0.0114308999                       2.481942
    ## Item3  0.0306677262                       2.991547
    ## Item4  0.0320145563                       3.482157
    ##       CategoricalVariableTimeOfDay_3 ContiniousVariable
    ## Item1                       3.030288          0.1997114
    ## Item2                       3.487059          0.3999603
    ## Item3                       3.963378          0.5994318
    ## Item4                       4.489226          0.7992137
    ## 
    ## 
    ## [[1]]$EMRepetitions
    ## [1] 6
    ## 
    ## [[1]]$last.loglik
    ## [1] -111428.9
    ## 
    ## [[1]]$nPara
    ## [1] 85
    ## 
    ## [[1]]$Sigma
    ## , , 1
    ## 
    ##           [,1]      [,2]      [,3]      [,4]
    ## [1,] 1.5043270 0.5095816 0.5275085 0.4806189
    ## [2,] 0.5095816 1.5078906 0.5317408 0.4927726
    ## [3,] 0.5275085 0.5317408 1.5341107 0.5158119
    ## [4,] 0.4806189 0.4927726 0.5158119 1.5036595
    ## 
    ## , , 2
    ## 
    ##           [,1]      [,2]      [,3]      [,4]
    ## [1,] 1.4764274 0.4891055 0.4956130 0.4745496
    ## [2,] 0.4891055 1.4769201 0.4882918 0.4775952
    ## [3,] 0.4956130 0.4882918 1.4914210 0.4712006
    ## [4,] 0.4745496 0.4775952 0.4712006 1.5006374
    ## 
    ## 
    ## [[1]]$LogLikelihood
    ##  [1] -141655.5 -133377.7 -118022.7 -111464.3 -111428.9 -111428.9        NA
    ##  [8]        NA        NA        NA        NA        NA        NA        NA
    ## [15]        NA        NA        NA        NA        NA        NA        NA
    ## [22]        NA        NA        NA        NA
    ## 
    ## [[1]]$EMiterationReset
    ##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [12] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [23] FALSE FALSE FALSE
    ## 
    ## [[1]]$PosteriorProbs
    ##        [,1] [,2]
    ##   [1,]    1    0
    ##   [2,]    1    0
    ##   [3,]    1    0
    ##   [4,]    0    1
    ##   [5,]    0    1
    ##   [6,]    0    1
    ##   [7,]    1    0
    ##   [8,]    1    0
    ##   [9,]    1    0
    ##  [10,]    1    0
    ##  [11,]    0    1
    ##  [12,]    1    0
    ##  [13,]    0    1
    ##  [14,]    1    0
    ##  [15,]    1    0
    ##  [16,]    0    1
    ##  [17,]    1    0
    ##  [18,]    1    0
    ##  [19,]    0    1
    ##  [20,]    0    1
    ##  [21,]    0    1
    ##  [22,]    1    0
    ##  [23,]    0    1
    ##  [24,]    1    0
    ##  [25,]    1    0
    ##  [26,]    1    0
    ##  [27,]    0    1
    ##  [28,]    1    0
    ##  [29,]    0    1
    ##  [30,]    0    1
    ##  [31,]    1    0
    ##  [32,]    0    1
    ##  [33,]    1    0
    ##  [34,]    0    1
    ##  [35,]    0    1
    ##  [36,]    1    0
    ##  [37,]    1    0
    ##  [38,]    0    1
    ##  [39,]    0    1
    ##  [40,]    0    1
    ##  [41,]    1    0
    ##  [42,]    0    1
    ##  [43,]    1    0
    ##  [44,]    0    1
    ##  [45,]    1    0
    ##  [46,]    1    0
    ##  [47,]    0    1
    ##  [48,]    0    1
    ##  [49,]    1    0
    ##  [50,]    1    0
    ##  [51,]    1    0
    ##  [52,]    1    0
    ##  [53,]    0    1
    ##  [54,]    1    0
    ##  [55,]    1    0
    ##  [56,]    0    1
    ##  [57,]    0    1
    ##  [58,]    1    0
    ##  [59,]    0    1
    ##  [60,]    0    1
    ##  [61,]    0    1
    ##  [62,]    1    0
    ##  [63,]    0    1
    ##  [64,]    1    0
    ##  [65,]    1    0
    ##  [66,]    0    1
    ##  [67,]    0    1
    ##  [68,]    1    0
    ##  [69,]    1    0
    ##  [70,]    1    0
    ##  [71,]    1    0
    ##  [72,]    0    1
    ##  [73,]    1    0
    ##  [74,]    1    0
    ##  [75,]    1    0
    ##  [76,]    1    0
    ##  [77,]    1    0
    ##  [78,]    1    0
    ##  [79,]    0    1
    ##  [80,]    1    0
    ##  [81,]    0    1
    ##  [82,]    1    0
    ##  [83,]    0    1
    ##  [84,]    0    1
    ##  [85,]    1    0
    ##  [86,]    0    1
    ##  [87,]    0    1
    ##  [88,]    0    1
    ##  [89,]    0    1
    ##  [90,]    1    0
    ##  [91,]    0    1
    ##  [92,]    1    0
    ##  [93,]    1    0
    ##  [94,]    0    1
    ##  [95,]    0    1
    ##  [96,]    0    1
    ##  [97,]    1    0
    ##  [98,]    0    1
    ##  [99,]    0    1
    ## [100,]    1    0
    ## [101,]    1    0
    ## [102,]    0    1
    ## [103,]    0    1
    ## [104,]    0    1
    ## [105,]    1    0
    ## [106,]    0    1
    ## [107,]    1    0
    ## [108,]    0    1
    ## [109,]    0    1
    ## [110,]    0    1
    ## [111,]    0    1
    ## [112,]    0    1
    ## [113,]    1    0
    ## [114,]    0    1
    ## [115,]    1    0
    ## [116,]    0    1
    ## [117,]    1    0
    ## [118,]    0    1
    ## [119,]    1    0
    ## [120,]    0    1
    ## 
    ## [[1]]$Lags
    ## [1] 1 1
    ## 
    ## [[1]]$Classification
    ##   [1] 1 1 1 2 2 2 1 1 1 1 2 1 2 1 1 2 1 1 2 2 2 1 2 1 1 1 2 1 2 2 1 2 1 2 2
    ##  [36] 1 1 2 2 2 1 2 1 2 1 1 2 2 1 1 1 1 2 1 1 2 2 1 2 2 2 1 2 1 1 2 2 1 1 1
    ##  [71] 1 2 1 1 1 1 1 1 2 1 2 1 2 2 1 2 2 2 2 1 2 1 1 2 2 2 1 2 2 1 1 2 2 2 1
    ## [106] 2 1 2 2 2 2 2 1 2 1 2 1 2 1 2
    ## 
    ## [[1]]$IC
    ## [1] 1.111175
    ## 
    ## [[1]]$SC
    ## [1] 1.135838
    ## 
    ## [[1]]$Proportions
    ##      .data_1 .data_2
    ## [1,]     0.5     0.5

``` r
Result$AllSolutions[[1]][[2]][[4]]
```

    ## $Converged
    ## [1] TRUE
    ## 
    ## $A
    ## , , 1
    ## 
    ##            Item1      Item2        Item3       Item4        Item1
    ## Item1  0.0106515  0.2671153  0.012149280 0.004481859 -0.018124224
    ## Item2 -0.4146793  0.3869261  0.264064519 0.300979921 -0.005133318
    ## Item3  0.6034681 -0.2958824 -0.006513921 0.280847617  0.024730087
    ## Item4  0.3108251 -0.3988043  0.392106037 0.554954640  0.007354099
    ##              Item2        Item3        Item4
    ## Item1  0.005448890  0.008323587 -0.008064031
    ## Item2  0.011581148 -0.004005090  0.004490731
    ## Item3 -0.018431289 -0.021319615  0.004455396
    ## Item4 -0.001376256  0.008704764 -0.015556552
    ## 
    ## , , 2
    ## 
    ##             Item1      Item2     Item3     Item4 Item1 Item2 Item3 Item4
    ## Item1  0.08479586 -0.1033597 0.1068198 0.2466591    NA    NA    NA    NA
    ## Item2 -0.21764302  0.1776732 0.1040631 0.1632604    NA    NA    NA    NA
    ## Item3  0.29650512 -0.3046279 0.1888623 0.1450859    NA    NA    NA    NA
    ## Item4  0.40965841 -0.4093675 0.1878278 0.4930011    NA    NA    NA    NA
    ## 
    ## 
    ## $B
    ## , , 1
    ## 
    ##          Intercept CategoricalVariableTimeOfDay_2
    ## Item1  0.015228873                       1.982025
    ## Item2  0.018958818                       2.452789
    ## Item3 -0.036192189                       3.013027
    ## Item4 -0.002321288                       3.494736
    ##       CategoricalVariableTimeOfDay_3 ContiniousVariable
    ## Item1                       2.967232          0.1000219
    ## Item2                       3.457505          0.2003108
    ## Item3                       4.033086          0.3009098
    ## Item4                       4.508692          0.3998918
    ## 
    ## , , 2
    ## 
    ##           Intercept CategoricalVariableTimeOfDay_2
    ## Item1 -0.0007224385                       2.001126
    ## Item2  0.0114308999                       2.481942
    ## Item3  0.0306677262                       2.991547
    ## Item4  0.0320145563                       3.482157
    ##       CategoricalVariableTimeOfDay_3 ContiniousVariable
    ## Item1                       3.030288          0.1997114
    ## Item2                       3.487059          0.3999603
    ## Item3                       3.963378          0.5994318
    ## Item4                       4.489226          0.7992137
    ## 
    ## 
    ## $EMRepetitions
    ## [1] 2
    ## 
    ## $last.loglik
    ## [1] -111067.5
    ## 
    ## $nPara
    ## [1] 101
    ## 
    ## $Sigma
    ## , , 1
    ## 
    ##           [,1]      [,2]      [,3]      [,4]
    ## [1,] 1.5055963 0.5099283 0.5263043 0.4817252
    ## [2,] 0.5099283 1.5101683 0.5321159 0.4944527
    ## [3,] 0.5263043 0.5321159 1.5338779 0.5166816
    ## [4,] 0.4817252 0.4944527 0.5166816 1.5049912
    ## 
    ## , , 2
    ## 
    ##           [,1]      [,2]      [,3]      [,4]
    ## [1,] 1.4764274 0.4891055 0.4956130 0.4745496
    ## [2,] 0.4891055 1.4769201 0.4882918 0.4775952
    ## [3,] 0.4956130 0.4882918 1.4914210 0.4712006
    ## [4,] 0.4745496 0.4775952 0.4712006 1.5006374
    ## 
    ## 
    ## $LogLikelihood
    ##  [1] -111067.5 -111067.5        NA        NA        NA        NA        NA
    ##  [8]        NA        NA        NA        NA        NA        NA        NA
    ## [15]        NA        NA        NA        NA        NA        NA        NA
    ## [22]        NA        NA        NA        NA
    ## 
    ## $EMiterationReset
    ##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [12] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [23] FALSE FALSE FALSE
    ## 
    ## $PosteriorProbs
    ##        [,1] [,2]
    ##   [1,]    1    0
    ##   [2,]    1    0
    ##   [3,]    1    0
    ##   [4,]    0    1
    ##   [5,]    0    1
    ##   [6,]    0    1
    ##   [7,]    1    0
    ##   [8,]    1    0
    ##   [9,]    1    0
    ##  [10,]    1    0
    ##  [11,]    0    1
    ##  [12,]    1    0
    ##  [13,]    0    1
    ##  [14,]    1    0
    ##  [15,]    1    0
    ##  [16,]    0    1
    ##  [17,]    1    0
    ##  [18,]    1    0
    ##  [19,]    0    1
    ##  [20,]    0    1
    ##  [21,]    0    1
    ##  [22,]    1    0
    ##  [23,]    0    1
    ##  [24,]    1    0
    ##  [25,]    1    0
    ##  [26,]    1    0
    ##  [27,]    0    1
    ##  [28,]    1    0
    ##  [29,]    0    1
    ##  [30,]    0    1
    ##  [31,]    1    0
    ##  [32,]    0    1
    ##  [33,]    1    0
    ##  [34,]    0    1
    ##  [35,]    0    1
    ##  [36,]    1    0
    ##  [37,]    1    0
    ##  [38,]    0    1
    ##  [39,]    0    1
    ##  [40,]    0    1
    ##  [41,]    1    0
    ##  [42,]    0    1
    ##  [43,]    1    0
    ##  [44,]    0    1
    ##  [45,]    1    0
    ##  [46,]    1    0
    ##  [47,]    0    1
    ##  [48,]    0    1
    ##  [49,]    1    0
    ##  [50,]    1    0
    ##  [51,]    1    0
    ##  [52,]    1    0
    ##  [53,]    0    1
    ##  [54,]    1    0
    ##  [55,]    1    0
    ##  [56,]    0    1
    ##  [57,]    0    1
    ##  [58,]    1    0
    ##  [59,]    0    1
    ##  [60,]    0    1
    ##  [61,]    0    1
    ##  [62,]    1    0
    ##  [63,]    0    1
    ##  [64,]    1    0
    ##  [65,]    1    0
    ##  [66,]    0    1
    ##  [67,]    0    1
    ##  [68,]    1    0
    ##  [69,]    1    0
    ##  [70,]    1    0
    ##  [71,]    1    0
    ##  [72,]    0    1
    ##  [73,]    1    0
    ##  [74,]    1    0
    ##  [75,]    1    0
    ##  [76,]    1    0
    ##  [77,]    1    0
    ##  [78,]    1    0
    ##  [79,]    0    1
    ##  [80,]    1    0
    ##  [81,]    0    1
    ##  [82,]    1    0
    ##  [83,]    0    1
    ##  [84,]    0    1
    ##  [85,]    1    0
    ##  [86,]    0    1
    ##  [87,]    0    1
    ##  [88,]    0    1
    ##  [89,]    0    1
    ##  [90,]    1    0
    ##  [91,]    0    1
    ##  [92,]    1    0
    ##  [93,]    1    0
    ##  [94,]    0    1
    ##  [95,]    0    1
    ##  [96,]    0    1
    ##  [97,]    1    0
    ##  [98,]    0    1
    ##  [99,]    0    1
    ## [100,]    1    0
    ## [101,]    1    0
    ## [102,]    0    1
    ## [103,]    0    1
    ## [104,]    0    1
    ## [105,]    1    0
    ## [106,]    0    1
    ## [107,]    1    0
    ## [108,]    0    1
    ## [109,]    0    1
    ## [110,]    0    1
    ## [111,]    0    1
    ## [112,]    0    1
    ## [113,]    1    0
    ## [114,]    0    1
    ## [115,]    1    0
    ## [116,]    0    1
    ## [117,]    1    0
    ## [118,]    0    1
    ## [119,]    1    0
    ## [120,]    0    1
    ## 
    ## $Lags
    ## [1] 2 1
    ## 
    ## $Classification
    ##   [1] 1 1 1 2 2 2 1 1 1 1 2 1 2 1 1 2 1 1 2 2 2 1 2 1 1 1 2 1 2 2 1 2 1 2 2
    ##  [36] 1 1 2 2 2 1 2 1 2 1 1 2 2 1 1 1 1 2 1 1 2 2 1 2 2 2 1 2 1 1 2 2 1 1 1
    ##  [71] 1 2 1 1 1 1 1 1 2 1 2 1 2 2 1 2 2 2 2 1 2 1 1 2 2 2 1 2 2 1 1 2 2 2 1
    ## [106] 2 1 2 2 2 2 2 1 2 1 2 1 2 1 2
    ## 
    ## $IC
    ## [1] 1.116614
    ## 
    ## $SC
    ## [1] 1.153753
    ## 
    ## $Proportions
    ##      .data_1 .data_2
    ## [1,]     0.5     0.5
