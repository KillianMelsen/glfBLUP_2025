###############################################################################
# Function that calculates the training set focal trait gBLUP based on
# multivariate data using a fast algorithm. Returns a vector of length n_train.
# Returns a list containing a matrix with training set gBLUPs for all traits and
# a named vector containing training set gBLUPs for just the focal trait.
# Arguments:
# Y           = n_train x p matrix of trait values
# K           = n x n kinship matrix (including training as well as test genotypes)
# X           = covariate design matrix, default is an n_train x 1 matrix with 1's for
#               the intercept.
# B           = covariate coefficient matrix, default is a 1 x p matrix with 0's 
#               (no intercepts)
# Vg          = p x p genetic covariance matrix
# Ve          = p x p residual covariance matrix
# C           = p x p inverse genetic covariance matrix
# D           = p x p inverse residual covariance matrix
# test.set    = vector of test genotype names
# train.set   = vector of training genotype names
# targetName  = name of the focal/target trait
###############################################################################

fast_o_gBLUP <- function(Y, K,
                   X = matrix(rep(1, nrow(K[train.set, train.set]))),
                   B = t(matrix(rep(0, ncol(Y)))),
                   Sg, Se, train.set, targetName) {
  
  p <- ncol(Y)
  Koo <- K[train.set, train.set]
  C <- solve(Sg)
  D <- MASS::ginv(Se) # Was solve()
  nc <- ncol(X)
  
  R <- solve(Koo)
  w <- eigen(R)
  Lambda.R <- diag(w$values)
  U <- w$vectors
  
  D.sqrt.inv <- MatrixRoot(MASS::ginv(D)) # Was solve()
  w <- eigen(D.sqrt.inv %*% C %*% D.sqrt.inv)
  Q.1 <- w$vectors
  Lambda.1 <- w$values
  
  if (nc > 0) {
    S.1     <- vec.inv.diag(x = Lambda.1, y = diag(Lambda.R)) * (t(U) %*% (Y - X %*% B) %*% MatrixRoot(D) %*% Q.1)
  } else {
    S.1     <- vec.inv.diag(x = Lambda.1, y = diag(Lambda.R)) * (t(U) %*% Y %*% MatrixRoot(D) %*% Q.1)
  }
  
  # D.sqrt.inv <- MatrixRoot(MASS::ginv(D)) # Was solve()
  # w <- eigen(D.sqrt.inv %*% C %*% D.sqrt.inv)
  # Q.1 <- w$vectors
  
  # All training gBLUPs:
  MU <- matrix(kronecker(D.sqrt.inv %*% Q.1, U) %*% matrix(S.1), ncol = p)
  colnames(MU) <- colnames(Y)
  rownames(MU) <- rownames(Y)
  focal <- MU[, targetName]
  focal <- as.numeric(focal)
  names(focal) <- rownames(MU)
  return(list(all_gBLUPs = MU,
              focal_gBLUP = focal))
}