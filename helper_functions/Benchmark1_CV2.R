Benchmark1_CV2 <- function(d, Sg, Se, K, r) {
  train.set <- as.character(unique(d$G[!is.na(d[,ncol(d)])]))
  test.set  <- as.character(setdiff(unique(d$G), train.set))
  targetName <- names(d)[length(names(d))]
  factors <- setdiff(colnames(Sg), targetName)
  # dm <- gfBLUPold:::genoMeans(data = d)
  # dm <- dm[match(unique(d$G), dm$G),]
  # rownames(dm) <- dm$G

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





