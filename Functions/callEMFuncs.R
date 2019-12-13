callEMFuncs <- function(Clusters, HighestLag, LowestLag, Rand, Rational, Initialization,
                        PreviousSol, IDNames, K, N, Y, X, Tni, qqq, nDepVar, 
                        PersStart, PersPDiffStart, PersEnd, PersStartU, PersEndU, Covariates,
                        Conv, it, val.init, ICType, smallestClN, SigmaIncrease)
{
    ### Loop over different K (# of clusters) values  -------------
    OutputAllK = vector("list", length(Clusters)) # with length(Clusters) many elements, containing the best solution for every K
    WHOLEOutputAllK  = vector("list", length(Clusters)) # with lenght(Clusters) many elements containing
    # the OutputListAllLags of every K, which contains all solutions for all Lags and all starts for that K
    
    ClustCount = 1
    for (K in Clusters)
    {
        cat(c("\n", K, "Clusters: "))
        cat(c("\n", K, "Clusters: "), file = "EMwarnings.txt", append = TRUE)
        
        LagsList = calculateLagList(K = K, HighestLag = HighestLag, LowestLag = LowestLag)
        LagCombinations = dim(LagsList)[1]
        
        # OutputListAllLags[[Lags]][[Start]]
        OutputListAllLags = createOutputList(LagCombinations = LagCombinations,
                                             Rand = Rand, Rational = Rational,
                                             Initialization = Initialization,
                                             PreviousSol = PreviousSol)
        # Fit[Lags, Start]
        FitAllLags = array(NA, dim = c(LagCombinations, Rand + as.numeric(Rational) + 
                                           as.numeric(!is.null(Initialization))
                                       + as.numeric(PreviousSol))) # to store fit of output
        
        CoeffsForRandoAndRationalList = callCalculateCoefficientsForRandoAndRational(
            Covariates = Covariates, K = K, N = N,
            nDepVar = nDepVar, qqq = qqq,
            HighestLag = HighestLag, LowestLag = LowestLag, 
            PersEnd = PersEnd, PersStart = PersStart, Y = Y,
            X = X, PersPDiffStart = PersPDiffStart)
        
        EMCallVec = c(as.character(1:Rand),
                      ifelse(Rational, "Rational", NULL),
                      switch(!is.null(Initialization), "FALSE" = "Initialization"), # default is NULL
                      ifelse(PreviousSol, "Previous", NULL))
        
        usePrevLagSol = FALSE   # Make sure the previous lag solution-using start is not called before a previous solution exists
        for (LagCounter in 1:LagCombinations)
        {
            
            cat(c('\n', " Lags:", LagsList[LagCounter, ], '\n'))
            cat(c('\n', " Lags:", LagsList[LagCounter, ], '\n'), file = "EMwarnings.txt", append = TRUE)
            ### Initialization Prerequesites: calcuate coefficeints passed to initial clustering solutions ###---------------------------------
            #PersPDiffStart, PersStartU etc can (must) be integers instead of vectors in calculateCoefficientsForRandoAndRational
            #### Random starts ###
            StartCounter = 0
            while (StartCounter != length(EMCallVec)) 
                # ToDo: in all EMFunc IDNames is passed but not used
            {
                
                StartCounter = StartCounter + 1
                cat(c("*", EMCallVec[StartCounter], "  "))
                cat(c("*", EMCallVec[StartCounter], "  "), file = "EMwarnings.txt", append = TRUE)
                OutputListAllLags[[LagCounter]][[StartCounter]] =
                    EMFunc(Init = EMInit(InitMT = 
                                             switch(EMCallVec[StartCounter],
                                                    "Rational" = InitRat(K = K,
                                                                         CoefficientsForRandoAndRational =
                                                                             CoeffsForRandoAndRationalList[[max(LagsList[LagCounter, ])]]),
                                                    
                                                    "Initialization" = val.init,
                                                    
                                                    "Previous" = if(usePrevLagSol)
                                                                {   ## Previous Sol ##
                                                                    t(dummy_cols(OutputListAllLags
                                                                     [[PrevBestRun[1]]][[PrevBestRun[2]]]$Classification,
                                                                     remove_first_dummy = FALSE)[ , -1])
                                                                }else
                                                                {   ## PseudoRand (a previous sol does not exist yet) ##
                                                                    InitPseudoRand(N = N, K = K, smallestClN = smallestClN,
                                                                                   CoefficientsForRandoAndRational =
                                                                                   CoeffsForRandoAndRationalList[[max(LagsList[LagCounter, ])]])
                                                                },
                                                    ### Default: PseudoRand initialization ###
                                                    InitPseudoRand(N = N, K = K, smallestClN = smallestClN,
                                                                   CoefficientsForRandoAndRational =
                                                                   CoeffsForRandoAndRationalList[[max(LagsList[LagCounter, ])]])
                                                    ),
                                         
                                         Y = Y, X = X, Lags = LagsList[LagCounter, ], K = K, N = N, Tni = Tni, qqq = qqq, nDepVar = nDepVar, 
                                         PersStart = PersStart, PersPDiffStart = PersPDiffStart, PersEnd = PersEnd,
                                         PersStartU = PersStartU, PersEndU = PersEndU, 
                                         Covariates = Covariates, smallestClN = smallestClN, SigmaIncrease = SigmaIncrease), 
                           IDNames = IDNames, Y = Y, X = X, K = K, N = N, Tni = Tni, qqq = qqq, nDepVar = nDepVar, 
                           PersPDiffStart = PersPDiffStart, PersEnd = PersEnd, PersStartU = PersStartU, PersEndU = PersEndU,
                           Covariates = Covariates, Conv = Conv, it = it, smallestClN = smallestClN, ICType = ICType, 
                           SigmaIncrease = SigmaIncrease)
                
                FitAllLags[LagCounter, StartCounter] = OutputListAllLags[[LagCounter]][[StartCounter]]$IC
               
            } # End of Start loop
            
            PrevBestRun = arrayInd(which.min(FitAllLags), dim(FitAllLags))
            usePrevLagSol = TRUE
            LagCounter = LagCounter + 1
            
        } # End of Lag loop
        
        # ToDo: tempo = proc.time() - ptm
        BestRunOneK = arrayInd(which.min(FitAllLags), dim(FitAllLags))
        # ToDo: update Classification with ID
        WHOLEOutputAllK[[ClustCount]] = OutputListAllLags
        OutputAllK[[ClustCount]] = OutputListAllLags[[BestRunOneK[1]]][[BestRunOneK[2]]]
        ClustCount = ClustCount + 1
        
    } # End of K loop
    
    invisible(list(BestSolutionsPerCluster = OutputAllK, AllSolutions = WHOLEOutputAllK))
} 
