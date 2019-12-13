# EuclideDistanceLag1 = matrix(c(.01, .05, .009, .049), 2, 2, byrow = TRUE) # first is used for medium cluster differences, second for large cluster differences
# EuclideDistanceLag2 = matrix(c(NA, NA, .001, .001), 2, 2, byrow = TRUE)  # first is used for medium cluster differences, second for large cluster differences

# EuclideDistanceLag1[i , m] # Cols gives distances, rows if the distance is for lag 1 or lag 2 ()
# EuclideDistanceLag2[i , m]


#### Caculate btwClDistMat used in Data Generation ####
CoefMatrices <- array(NA, dim = c(numberVariables, numberVariables, max(nClusters)))


CoefMatrices[ , , 1] = matrix(rep(0, numberVariables * numberVariables),
                              ncol = numberVariables, nrow = numberVariables, byrow = TRUE)


CoefMatrices[ , , 2] = matrix(c(1, 0, 0, 1,
                                1, 0, 1, 0,
                                1, 0, -1, 0,
                                1, -1, 0, 0),
                              ncol = numberVariables, nrow = numberVariables, byrow = TRUE)

CoefMatrices[ , , 3] = matrix(c(1, 0, 1, 0,
                                1, 0, 0, 1,
                                1, -1, 0, 0,
                                1, 0, -1, 0),
                              ncol = numberVariables, nrow = numberVariables, byrow = TRUE)

CoefMatrices[ , , 4] = matrix(c(0, 0, 1, 1,
                                1, -1, 0, 0,
                                1, 0, 0, 1,
                                0, -1, -1, 0),
                              ncol = numberVariables, nrow = numberVariables, byrow = TRUE)

btwClDistMat <- array(NA, dim = c(numberVariables, max(nLags)*numberVariables, max(nClusters), length(nLags)))

### Lag 1 btwClDistMat[ , , , 1] ###
btwClDistMat[ , , 1, 1] = cbind(CoefMatrices[ , , 1], CoefMatrices[ , , 1])
# The first btwClDistMat is filled with zeros, is for first matrix which is not altered
btwClDistMat[ , , 2, 1] = cbind(CoefMatrices[ , , 2], CoefMatrices[ , , 1])
btwClDistMat[ , , 3, 1] = cbind(CoefMatrices[ , , 3], CoefMatrices[ , , 1])
btwClDistMat[ , , 4, 1] = cbind(CoefMatrices[ , , 4], CoefMatrices[ , , 1])

### Lag 2 btwClDistMat[ , , , 2]
btwClDistMat[ , , 1, 2] = cbind(CoefMatrices[ , , 1], CoefMatrices[ , , 1])
# The first btwClDistMat is filled with zeros, is for first matrix which is not altered
btwClDistMat[ , , 2, 2] = cbind(CoefMatrices[ , 1:2, 2], CoefMatrices[ , 1:2, 1], CoefMatrices[ , 3:4, 2], CoefMatrices[ , 3:4, 1])
btwClDistMat[ , , 3, 2] = cbind(CoefMatrices[ , 1:2, 3], CoefMatrices[ , 1:2, 1], CoefMatrices[ , 3:4, 3], CoefMatrices[ , 3:4, 1])
btwClDistMat[ , , 4, 2] = cbind(CoefMatrices[ , 1:2, 4], CoefMatrices[ , 1:2, 1], CoefMatrices[ , 3:4, 4], CoefMatrices[ , 3:4, 1])

