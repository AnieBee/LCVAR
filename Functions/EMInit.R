
EMInit <- function(InitMT, Y, X, Lags, K, N, Tni, qqq, nDepVar, PersStart, PersPDiffStart, PersEnd, PersStartU, PersEndU, Covariates, smallestClN, SigmaIncrease)
    # Calls determineLagOrder and reorderLags to switch Lags
{
    ListCCC = checkComponentsCollapsed(K = K, N = N, FZY = t(InitMT),
                                       smallestClN = smallestClN,
                                       EMiteration = 0, crisp = TRUE)
    memb = t(ListCCC$FZY)
    stopifnot( ! any(colSums(memb) != 1))
    # Do the EM Initialization given memb and tau as determined by either InitRat or InitRand or given by val.init
    ### initialization of weights (priors) based on assignment
    tau = t(rowSums(memb) / N) 
    # tau is not needed in the EMInit func but is passed to EMfunc via retrun statement
    
    # Constraints on B
    DimensionsBasedonConstraints = constraintsOnB(Covariates, K, N)
    
    # Partition matrices
    A = array(NA, dim = c(nDepVar, nDepVar * max(Lags), K))
    UZero = array(NA, dim = c(nDepVar, PersEndU[[min(Lags)]][N], K)) 
    # Vector of u_{ikt}s # U is not of same length as Y, Y contains N many Lags*m pre-samles
    # U has different lenght for different clusters when clusters have different lag numbers
    Sigma = array(0, dim = c(nDepVar, nDepVar, K)) 
    B = array(0, dim = c(nDepVar, qqq, DimensionsBasedonConstraints$BNumbVersions)) 
    # B is of dim(3) = K if B are unequal per cluster, if B are equal is of dim(3) = 1, dim(3) = N if individual-specific
    
    ### Calculate B estimate using LS estimation & Calculate WZero = W(0) based on Ls estimated B -----------------
    BZero = array(0, dim = c(nDepVar, qqq, DimensionsBasedonConstraints$BNumbVersions)) 
    # B is of dim(3) = K if B are unequal per cluster, if B are equal is of dim(3) = 1, dim(3) = N if individual-specific
    WZero = array(0, dim = c(nDepVar, PersEnd[N], 1)) # WZero is two dimensional NO MATTER the B constraint (Wk in EM is 3 dimensional) 
    # because every person belongs to only a single cluster, for every person there is only a single W calculated
    
    WZero = calculateBandWZero(Covariates = Covariates, K = K, memb = memb, Y = Y, X = X ,
                               PersStart = PersStart, PersEnd = PersEnd, BZero = BZero, WZero = WZero)
    
    ###### Reorder Lags depending on order exhibited in ClusterVARcoeffs -----------------
    ClusterVARcoeffs = calculateA(K = K, WkNumbVersions = 1, 
                                  # One Wzero exists for all clusters, because memb is 0 or 1 set WkNumbVersions
                                  # to 1 because WZero instead of Wk is passed here
                                  N = N, Wk = WZero, PersPDiffStart = PersPDiffStart, 
                                  PersEnd = PersEnd, Lags = rep(max(Lags), K), FZY = t(memb),
                                  A = A, nDepVar)
    
    Lags = reorderLags(Lags = Lags, K = K,
                         newOrder = determineLagOrder(Lags = Lags, K = K,
                                                      ClusterVARcoeffs = ClusterVARcoeffs,
                                                      nDepVar = nDepVar))
    ###### Calculate A (based on Lags), UZero, S and B based on W(0) = WZero --------------
    A = calculateA(K = K, WkNumbVersions = 1,
                   # One Wzero exists for all clusters, set WkNumbVersions to 1 because WZero instead of Wk is passed here
                   N = N, Wk = WZero, PersPDiffStart = PersPDiffStart, 
                   PersEnd = PersEnd, Lags = Lags, FZY = t(memb), A = A,
                   nDepVar = nDepVar)
    
    UZero = calculateU(K = K, WkNumbVersions = 1, N = N, PersPDiffStart = PersPDiffStart,
                       PersEnd = PersEnd, U = UZero, Wk = WZero, A = A, Lags = Lags,
                       nDepVar = nDepVar)
    
    Sigma = calculateSigma(K = K, N = N, FZY = t(memb), U = UZero, PersStartU = PersStartU,
                           PersEndU = PersEndU, Tni = Tni, Sigma = Sigma, Lags = Lags)
    SigmaList = checkSingularitySigma(nDepVar = nDepVar, K = K, Sigma = Sigma, EMiteration = 0)
    Sigma = SigmaList$Sigma
    Sigma[ , , ListCCC$resetCl] = Sigma[ , , ListCCC$resetCl] + SigmaIncrease
    
    ## Calculate B(0) depending on Covariate constraint----------------------------
    # B(0) calculation differs from B cacluation in using memb[j, i] instead of f.z.y[ i, j]
    B = calculateB(Covariates = Covariates, K = K, nDepVar = nDepVar, A = A,
                   Sigma = Sigma, N = N, PersPDiffStart = PersPDiffStart, 
                   PersEnd = PersEnd, X = X, Y = Y, Lags = Lags,
                   FZY = t(memb), qqq = qqq, B = B)
    
    #### Return initialization to EM call
    invisible(list(A = A, B = B, Sigma = Sigma, tau = tau, Lags = Lags))
    
}
