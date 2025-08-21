# This file contains a number of helper functions that are sourced by some of the
# other scripts. Each function, its arguments, and returned value is described
# in a comment above the function.


# Benchmark1_CV2 produces benchmark predictions for simulated data
# given the true, simulated covariance matrices in the CV2 setting.
#
# d: A dataframe with the genotype factor in the first column, the focal trait
#    in the final column (with values for the test set as NA), and factor scores inbetween.
# Sg: The true, simulated genetic covariance matrix of factors and focal trait.
# Se: The true, simulated residual covariance matrix of factors and focal trait.
# K: The kinship matrix.
# r: The number of replicates.
#
# This function returns a list with benchmark predictions.
Benchmark1_CV2 <- function(d, Sg, Se, K, r) {
  train.set <- as.character(unique(d$G[!is.na(d[,ncol(d)])]))
  test.set  <- as.character(setdiff(unique(d$G), train.set))
  targetName <- names(d)[length(names(d))]
  factors <- setdiff(colnames(Sg), targetName)
  
  BLUPs_train <- fast_o_gBLUP(Y = as.matrix(d[train.set, -1]),
                              K = K, Sg = Sg, Se = Se/r, train.set = train.set,
                              targetName = targetName)
  
  BLUPs_test <- train2test(BLUPs_train$all_gBLUPs, K = K, train.set = train.set,
                           test.set = test.set)
  
  BLUP_Y_test <- BLUPs_test[, targetName]
  BLUP_S_test <- BLUPs_test[, setdiff(colnames(BLUPs_test), targetName)]
  Ktt <- K[test.set, test.set]
  Kto <- K[test.set, train.set]
  Koo <- K[train.set, train.set]
  
  Vc <- kronecker(Sg[factors, factors], solve(Ktt)) +
    kronecker(Se[factors, factors] / r, diag(dim(Ktt)[1]))
  
  focalBLUP_test <- as.numeric(matrix(BLUP_Y_test) +
                                 kronecker(t(as.matrix(Sg[factors, targetName])), solve(Ktt)) %*%
                                 solve(Vc) %*%
                                 (matrix(as.matrix(d[test.set, factors])) - matrix(BLUP_S_test)))
  
  names(focalBLUP_test) <- test.set
  return(list(preds = focalBLUP_test))
}

# Benchmark1 produces benchmark predictions for simulated data
# given the true, simulated covariance matrices in the CV1 setting.
#
# d: A dataframe with the genotype factor in the first column, the focal trait
#    in the final column (with values for the test set as NA), and factor scores inbetween
#    (factor scores are also set to NA for the test set in CV1).
# Sg: The true, simulated genetic covariance matrix of factors and focal trait.
# Se: The true, simulated residual covariance matrix of factors and focal trait.
# K: The kinship matrix.
# r: The number of replicates.
#
# This function returns a list with benchmark predictions.
benchmark1 <- function(d, Sg, Se, K, r) {
  train.set <- as.character(unique(d$G[!is.na(d[,ncol(d)])]))
  test.set  <- as.character(setdiff(unique(d$G), train.set))
  targetName <- names(d)[length(names(d))]
  d.train <- d[which(d$G %in% train.set), ]
  
  focalBLUP_train <- fast_o_gBLUP(Y = as.matrix(d.train[, -1]), K = K, Sg = Sg,
                                  Se = Se / r, train.set = train.set,
                                  targetName = targetName)
  
  focalBLUP_test <- train2test(focalBLUP_train$focal_gBLUP, K = K,
                               train.set = train.set, test.set = test.set)
  
  return(list(preds = focalBLUP_test))
}

# estimate_gcor takes observed trait values and genomic predictions and estimates
# the genetic correlation between them. This fulfills the parametric method for
# evaluation genomic prediction accuracy from Runcie and Cheng 2019.
#
# data: data.frame with columns 'ID','obs','pred'.
# Knn: kinship among the testing set, all data$ID must be in rownames(Knn).
# sKnn: optional; result of \code{svd(Knn)}.
# method: analysis method to use.
# normalize: Should \code{obs} and \code{pred} be standardized to unit variance?
# control: List of optional arguments to modify behavior of methods.
#
# This function returns a named vector with the estimated genetic correlation and
# (if using method = "MCMCglmm"), the R-hat MCMC convergence diagnostic.
estimate_gcor = function(data, Knn, sKnn = NULL, method = c('MCMCglmm','sommer'), normalize = c(T,F), control = c()) {
  if(normalize) {
    data$obs = scale(data$obs)
    data$pred = scale(data$pred)
  }
  
  result = switch(method,
                  MCMCglmm = {
                    require(MCMCglmm)
                    if(is.null(sKnn)) sKnn = svd(Knn)
                    
                    Dinv = as(diag(1/sKnn$d),'dgCMatrix')
                    rownames(Dinv) = colnames(Dinv) = rownames(Knn)
                    fixed = formula('cbind(ut_obs,ut_pred) ~ 0+ut1:trait')
                    
                    n_test = nrow(data)
                    data$ut1 = t(sKnn$u) %*% matrix(1,nr = n_test,nc=1)
                    data$ut_obs = t(sKnn$u) %*% data$obs
                    data$ut_pred = t(sKnn$u) %*% data$pred
                    
                    prior = list(R = list(V = diag(c(.5,.01),2),nu = 3),G = list(G1 = list(V = diag(c(.5,.5),2),nu = 3,alpha.mu = rep(0,2),alpha.V = diag(1,2))))#
                    nitt = 30000
                    thin = 50
                    burn = 5000
                    if('prior' %in% names(control)) prior = control$prior
                    if('nitt' %in% names(control)) nitt = control$nitt
                    if('thin' %in% names(control)) thin = control$thin
                    if('burn' %in% names(control)) burn = control$burn
                    
                    m_gcor <- MCMCglmm(fixed,random = ~us(trait):ID,rcov = ~us(trait):units,data = data,prior = prior,family = rep('gaussian',2),pl=F,pr=T,
                                       ginverse = list(ID = Dinv),nitt = nitt,burn = burn,thin = thin,verbose = F)
                    
                    h2_obs = apply(m_gcor$VCV,1,function(x) (x[1]/(x[1] + x[5])))
                    h2_pred = apply(m_gcor$VCV,1,function(x) (x[4]/(x[4] + x[8])))
                    g_cor = apply(m_gcor$VCV,1,function(x) x[2]/sqrt(x[1]*x[4]))
                    r_cor = apply(m_gcor$VCV,1,function(x) x[6]/sqrt(x[5]*x[8]))
                    
                    return(c(g_cor = mean(g_cor * sqrt(h2_pred)), Rhat = rstan::Rhat(g_cor * sqrt(h2_pred))))
                    
                  },
                  sommer = {
                    require(sommer)
                    
                    data2 = rbind(data.frame(data,y = data$obs,Type = 'obs'),
                                  data.frame(data,y=data$pred,Type='pred'))
                    data2 = droplevels(data2)
                    
                    g_cor3 = try({
                      ms = mmer(y ~ 0 + Type,
                                random = ~ vsr(usr(Type), ID, Gu = Knn),
                                rcov = ~ vsr(usr(Type), units),
                                data = data2,
                                verbose = 0)
                      pin(ms,g_cor~V2/sqrt(V1*V3) * sqrt(V3/(V3+V6)))
                    },silent = T)
                    if(class(g_cor3) == 'try-error') {
                      ms = mmer(y ~ 0 + Type,
                                random = ~ vsr(usr(Type), ID, Gu = Knn),
                                rcov = ~ vsr(dsr(Type), units),
                                data = data2,
                                verbose = 0)
                      # V2 is genetic covariance between y and u, V1 and V3 are genetic variances of y and u, so
                      # V2/sqrt(V1*V3) = genetic correlation between y and u. V3 is genetic variance of u and V5
                      # is residual variance of u, so V3/(V3+V5) is heritability of u.
                      g_cor3 = try(pin(ms,g_cor~V2/sqrt(V1*V3) * sqrt(V3/(V3+V5))),silent = T)
                      if(class(g_cor3) == 'try-error') g_cor3 = c(Estimate = NA,SE = NA)
                    }
                    return(unlist(g_cor3))
                    
                  }
  )
  return(result)
}

# fast_o_gBLUP is a function that computes training set BLUPs based on multivariate
# data using a fast algorithm.
#
# Y: n_train x p matrix of trait values.
# K: n x n kinship matrix (including training as well as test genotypes).
# X: covariate design matrix, default is an n_train x 1 matrix with 1's for the intercept.
# B: covariate coefficient matrix, default is a 1 x p matrix with 0's (no intercepts).
# Sg: p x p genetic covariance matrix.
# Se: p x p residual covariance matrix.
# train.set: vector of training genotype names.
# targetName: name of the focal/target trait.
#
# This functions returns a list with a matrix containing BLUPs for all traits and
# a vector containing the focal trait BLUPs only.
fast_o_gBLUP <- function(Y, K,
                         X = matrix(rep(1, nrow(K[train.set, train.set]))),
                         B = t(matrix(rep(0, ncol(Y)))),
                         Sg, Se, train.set, targetName) {
  
  p <- ncol(Y)
  Koo <- K[train.set, train.set]
  C <- solve(Sg)
  D <- MASS::ginv(Se)
  nc <- ncol(X)
  
  R <- solve(Koo)
  w <- eigen(R)
  Lambda.R <- diag(w$values)
  U <- w$vectors
  
  D.sqrt.inv <- MatrixRoot(MASS::ginv(D))
  w <- eigen(D.sqrt.inv %*% C %*% D.sqrt.inv)
  Q.1 <- w$vectors
  Lambda.1 <- w$values
  
  if (nc > 0) {
    S.1     <- vec.inv.diag(x = Lambda.1, y = diag(Lambda.R)) * (t(U) %*% (Y - X %*% B) %*% MatrixRoot(D) %*% Q.1)
  } else {
    S.1     <- vec.inv.diag(x = Lambda.1, y = diag(Lambda.R)) * (t(U) %*% Y %*% MatrixRoot(D) %*% Q.1)
  }
  
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

# MatrixRoot is a function that computes the square roots of a matrix.
#
# x: The matrix of which the square roots must be computed.
#
# This function returns a matrix that is the square root of x.
MatrixRoot <- function(x) {
  x.eig <- eigen(as.matrix(x), symmetric = TRUE)
  if (length(x.eig$values) > 1) {
    x.sqrt <- x.eig$vectors %*% diag(sqrt(x.eig$values)) %*% solve(x.eig$vectors)
  } else {
    x.sqrt <- x.eig$vectors %*% matrix(sqrt(x.eig$values)) %*% solve(x.eig$vectors)
  }
  return(x.sqrt)
}

# train2test is a function that converts BLUPs for a training set to BLUPs for
# a test set using the kinship matrix.
#
# BLUPs: Vector or matrix of training set BLUPs to be used for calculating the
#        test set BLUPs.
# K: Kinship that contains the training and test genotypes.
# train.set: Vector containing training set genotype names.
# test.set: Vector containing test set genotype names.
#
# This function returns a vector or matrix of test set BLUPs.
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

# vec.inv.diag is a function that directly computes the trace of Lambda_1^{*} and
# Lambda_2^{*}.
#
# x: A numeric vector
# y: A numeric vector
#
# This function returns a matrix containing the traces of Lambda_1^{*} and
# Lambda_2^{*}.
vec.inv.diag <- function(x, y) {
  p <- length(x)
  n <- length(y)
  z <- sapply(x, function(x_i) {return(1 / (rep(1, n) + x_i * y))})
  return(z)
}



