
EMFunc <- function(Init, IDNames, Y, X, K, N, Tni, qqq, nDepVar, 
                   PersPDiffStart, PersEnd, PersStartU, PersEndU, Covariates,
                   Conv, it, smallestClN, ICType, SigmaIncrease)
{

  likelihood = rep(NA, it)
  iterationReset = rep(FALSE, it)         
  # Have any likelihoods/posteriors been reset in this iteration (FALSE means no reset in this iteration)
  EMiteration = 0 # counts the iteration
  ratio = 1000
  lik = -10000000
  lowest.Likelihood = -9000
  
  FYZ = matrix(0, N, K) # prob(W_{i}| Sigma_k, A_K, B_K), Normal weighted by prior 
  FZY = matrix(0, N, K) # posterior group membership (pi_{ik})
  A = Init$A
  B = Init$B
  Sigma = Init$Sigma # Sigma has to be checked for singularity, is checked for it at end of EMInit
  tau = Init$tau
  Lags = Init$Lags # contains switched Lags order to match lag order the clusters exhibited in EMInit
  U = array(NA, dim = c(nDepVar, PersEndU[[min(Lags)]][N], K))
  # Vector of u_{ikt}s # U is not of same length as Y, Y contains N many Lags*m pre-samples

  # Constraints on B #
  DimensionsBasedonConstraints = constraintsOnB(Covariates, K, N)
  Wk = array(0, dim = c(nDepVar, PersEnd[N], DimensionsBasedonConstraints$WkNumbVersions))
  # Wk is of dim(3) = K if B are unequal, if B are equal OR  individual-specific  dim(3) = 1
  
  while ( (EMiteration < it) & (ratio > Conv) )
  {#### E & M steps -------------
      
      EMiteration = EMiteration + 1

      ########  E-STEP -----------------------------------------
      
      ## Check: Sigma ----------------
      Wk = calculateW(Covariates = Covariates, K = K, Wk = Wk, Y = Y, B = B, X = X)
      
      U = calculateU(K = K, WkNumbVersions = DimensionsBasedonConstraints$WkNumbVersions,
                     N = N, PersPDiffStart = PersPDiffStart,
                     PersEnd = PersEnd, U = U, Wk = Wk, A = A,
                     Lags = Lags, nDepVar = nDepVar)
      
      # calculate the FYZs -----------
      # all FYZs and FZYs are in log form to avoid underflow
      # E step, or dvmnorm could be implemented in Rcpp
      FYZ = calculateFYZ(K = K, N = N, FYZ = FYZ, U = U, PersStartU = PersStartU,
                         PersEndU = PersEndU, nDepVar = nDepVar, Sigma = Sigma,
                         Lags = Lags)
      
      ## Check: likelihoods ----------
      FYZListCL = checkLikelihoodsNA(FYZ = FYZ, EMiteration = EMiteration)
      FYZ = FYZListCL$FYZ
      
      # calculate FZY
      calcPostList = calculatePosterior(N = N, K = K, tau = tau, FYZ = FYZ,
                                        lowest.Likelihood = lowest.Likelihood,
                                        EMiteration = EMiteration, FZY = FZY,
                                        UseFZY = TRUE)
      FZY = calcPostList$FZY
      # calculate FZY = posterior(pi_{ik})
      
      
      ### Checks --------------------
      ## Check for outliers
      FZYListCO = checkOutliers(K = K, FZY = FZY, EMiteration = EMiteration)
      
      ## Check posteriors are not NA ##
      FZYListCP = checkPosteriorsNA(FZY = FZYListCO$FZY, K = K, EMiteration = EMiteration)
      
      ## Check no component collapses onto a single point 
      FZYListCCC = checkComponentsCollapsed(K = K, N = N, FZY = FZYListCP$FZY,
                                            smallestClN = smallestClN,
                                            EMiteration = EMiteration)
      FZY = FZYListCCC$FZY
      
      ########  M-STEP    ########  
      ### calculate A -----------------
      A = calculateA(K = K, WkNumbVersions = DimensionsBasedonConstraints$WkNumbVersions,
                     N = N, Wk = Wk, PersPDiffStart = PersPDiffStart,
                     PersEnd = PersEnd, Lags = Lags, FZY = FZY,
                     A = A, nDepVar = nDepVar)
      # those places of A that will not be filled (because of lower lag number in some clusters)
      # will always be zero anyway because they are never filled, once A is calculated
      # in EMFunc, cannot be changed back
      
      ### calculate Sigma (S) ---------------------
      Sigma = calculateSigma(K = K, N = N, FZY = FZY, U = U, PersStartU = PersStartU,
                             PersEndU = PersEndU, Tni = Tni, Sigma = Sigma, Lags = Lags)
      SigmaList = checkSingularitySigma(nDepVar = nDepVar, K = K, Sigma = Sigma, EMiteration = EMiteration)
      Sigma = SigmaList$Sigma
      Sigma[ , , FZYListCCC$resetCl] = Sigma[ , , FZYListCCC$resetCl] + SigmaIncrease # Increase variance of components indicated by FZYListCCC 
      
      ## Calculate B depending on Covariate constraint -------------
      B = calculateB(Covariates = Covariates, K = K, nDepVar = nDepVar, A = A,
                     Sigma = Sigma, N = N, PersPDiffStart = PersPDiffStart,
                     PersEnd = PersEnd, X = X, Y = Y, Lags = Lags, FZY = FZY,
                     qqq = qqq, B = B)
      
      ### Calculate tau -------------------
      tau = calculateTau(K = K, tau = tau, FZY = FZY, N = N)
      
      ### end of M-step ###
      
      #### Calculate current log likelihood (temp) and compare to old likelihood (lik) ---
      # calculating the denominator of pi _{ik} using updated tau
      calcPostList2 = calculatePosterior(N = N, K = K, tau = tau, FYZ = FYZ,
                                         lowest.Likelihood = lowest.Likelihood,
                                         EMiteration = EMiteration, FZY = FZY,
                                         UseFZY = FALSE)
      likelihood[EMiteration] = calcPostList2$logLikelihood
      
      ### Has EM been reset in this iteration? --------------
      iterationReset[EMiteration] =  FYZListCL$iterationReset |
          calcPostList$iterationReset |
          FZYListCP$iterationReset | FZYListCCC$iterationReset |
          SigmaList$iterationReset | FZYListCO$iterationReset |
          calcPostList2$iterationReset 
      
      
      
      ## compare to old log likelihood (lik)
      ratio = calculateRatio(temp = likelihood[EMiteration], lik = lik,
                             EMiteration = EMiteration, Conv = Conv,
                             L2iterationReset = ifelse(EMiteration > 1,
                                                       iterationReset[EMiteration - 1] | iterationReset[EMiteration],
                                                       iterationReset[EMiteration]))
      lik = likelihood[EMiteration] 
  } ### end of while loop (end of E & M steps)-------
  
  Converged = checkConvergence(Conv = Conv, ratio = ratio,
                               L2iterationReset = iterationReset[EMiteration - 1]
                               | iterationReset[EMiteration])
  nPara = calculateNPara(Lags = Lags, nDepVar = nDepVar, K = K,
                         BnumbVersions = DimensionsBasedonConstraints$BNumbVersions,
                         ncovariates = qqq)
  IC = calculateIC(ICType = ICType, Sigma = Sigma, Lags = Lags,
                      nDepVar = nDepVar, K = K, N = N, FZY = FZY,
                      Tni = Tni, tau = tau)

  SC = calculateIC(ICType = "SC", Sigma = Sigma, Lags = Lags,
                      nDepVar = nDepVar, K = K, N = N, FZY = FZY,
                      Tni = Tni, tau = tau)
  
  Classification = apply(FZY, 1, which.max)
  last.lik = likelihood[EMiteration] 
  
  # Use ID names to return Classification and user knows what classification means
  
  ## the below is not implemented yet, Data is not passed to this function
  #Classification <- cbind(IDNames, Classification)
  #colnames(Classification) <- c("ID Name", "Cluster")
  
  colnames(B) = rownames(X) # Name every col in Array with Covariates Variable name
  rownames(B) = rownames(Y) # Name every row in B with Endogenous Variable Names
  colnames(A) = rep(rownames(Y), max(Lags)) 
  rownames(A) = rownames(Y)
  
  #Intercept is the "reference" group in case categorical variables are included
  invisible(list(Converged = Converged, A = A, B = B, EMRepetitions = EMiteration, 
                 last.loglik = last.lik, nPara = nPara, Sigma = Sigma, LogLikelihood = likelihood,
                 EMiterationReset = iterationReset, PosteriorProbs = FZY, Lags = Lags,
                 Classification = Classification,
                 IC = IC, SC = SC,  Proportions = tau))
  
} # end of EMfunc
