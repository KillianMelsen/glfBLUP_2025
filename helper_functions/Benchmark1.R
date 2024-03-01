benchmark1 <- function(d, Sg, Se, K, r) {
  train.set <- as.character(unique(d$G[!is.na(d[,ncol(d)])]))
  test.set  <- as.character(setdiff(unique(d$G), train.set))
  targetName <- names(d)[length(names(d))]
  # dm <- getMeans(suffStat = d)
  # dm <- dm[which(dm$G %in% train.set), ]
  d.train <- d[which(d$G %in% train.set), ]

  focalBLUP_train <- fast_o_gBLUP(Y = as.matrix(d.train[, -1]), K = K, Sg = Sg,
                                  Se = Se / r, train.set = train.set,
                                  targetName = targetName)

  focalBLUP_test <- train2test(focalBLUP_train$focal_gBLUP, K = K,
                               train.set = train.set, test.set = test.set)

  return(list(preds = focalBLUP_test))
}
