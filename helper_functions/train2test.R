###############################################################################
# Function that converts focal trait BLUPs for a training set to focal trait
# BLUPs for the test set. Returns a named vector of test set gBLUPs.
# Arguments:
# BLUPs     = Vector of training set BLUPs to be used for calculating the test
#             set BLUPs.
# K         = Kinship that contains the training and test genotypes.
# train.set = Vector containing training set genotype names.
# test.set  = Vector containing test set genotype names.
###############################################################################

train2test <- function(BLUPs, K, train.set, test.set) {
  Koo <- K[train.set, train.set]
  Kto <- K[test.set, train.set]
  focalBLUP_test <- Kto %*% solve(Koo) %*% BLUPs
  
  if (is.null(dim(BLUPs))) {
    focalBLUP_test <- as.numeric(focalBLUP_test)
    names(focalBLUP_test) <- test.set
  } else {
    rownames(focalBLUP_test) <- test.set
  }
  return(focalBLUP_test)
}



