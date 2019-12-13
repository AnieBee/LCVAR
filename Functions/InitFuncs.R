

### Rational ###
InitRat <- function(K, CoefficientsForRandoAndRational){
  #if (K > 1) memb <- kmeans(t(CoefficientsForRandoAndRational), K)$cl else memb <- rep(1, N)
    classification = kmeans(t(CoefficientsForRandoAndRational), K)$cl
    memb = t(dummy_cols(factor(classification, levels = as.character(1:K)), remove_first_dummy = FALSE)[ , -1])
  
  invisible(memb)
  
}

### Random Centres and Euclidean Distance ### 
# Random initialization used by Michael16: select randomly k cluster centres and partition other observations in that have shortest Euclidean distance
InitPseudoRand <- function(N, K, smallestClN, CoefficientsForRandoAndRational){
  # Random initialization: splitting people randomly in groups (memb gives random group membership)
  EuclideanDistanceFromCenters = array(NA, dim = c(N, K))
  memb = t(rep(1, N)) # initialize memb with 1s so the while loop
  # is entered because R does not have do-while statements
  
  while(sum(table(factor(apply(memb, 2, which.max), levels = as.character(1:K))) < smallestClN))
      # make sure no cluster contains less than smallestClN persons
  { 
    ClustCenters = sample(1:N, K, replace = FALSE) 
    # Randomly select K individuals as cluster centers
    for(i in 1:N)
    {
      for(j in 1:K)
      {
        EuclideanDistanceFromCenters[i, j] = dist(rbind(CoefficientsForRandoAndRational[, i],
                                                        as.vector(CoefficientsForRandoAndRational[ , ClustCenters[j]])),
                                                  method = "euclidean")
        # Above is equivalent to: EuclideanDistanceFromCenters[i, j] <- 
        # sqrt(sum((CoefficientsForRandoAndRational[ , i] - as.vector(CoefficientsForRandoAndRational[ , ClustCenters[j]]))^2))
      }
    }
    memb = t(dummy_cols(factor(apply(EuclideanDistanceFromCenters, 1, which.min),
                               levels = as.character(1:K)), remove_first_dummy = FALSE)[ , -1]) 
    # memb is of dimension K x N
    
  }
  
  invisible(memb)
  
}

### Random ###
InitRand <- function(N, K)
{
  
    #Random initialization: splitting people randomly in groups (memb gives random group membership)
    memb = t(rep(0, N)) # initialize memb with zero so the while loop is entered because R does not have do-while statements
    while(sum(table(factor(apply(memb, 2, which.max), levels = as.character(1:K))) < smallestClN))
        # make sure no cluster contains less than smallestClN persons
    {
        memb = rmultinom(N, size = 1, prob = rep(1 / K, K)) 
    }
  
    ### Alternative initialization of priors #might lead on average to more extreme priors, might be good
    # tau = runif(K, min = .1, max = 1)
    # tau = t(tau/sum(tau))
  
    invisible(memb)

}