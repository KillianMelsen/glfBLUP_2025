start <- Sys.time()

# Loading libraries:
library(rlist)
library(tictoc)

# Setting seed:
set.seed(1997)

# Total number of secondary features:
ps <- seq(100, 2500, 100)

# Loading kinship:
load("genotypes/K_sim.RData"); rm(M)

n.train <- 300

p <- ps[1]
for (p in ps) {

  tic(sprintf("p = %s", p))
  
  # Loading raw simulated dataset:
  datalist <- list.load(file = sprintf("timing/datasets/timing_p%d.RData", p))

  # Randomly choosing training and test genotypes:
  train.set <- sample(rownames(K), n.train)
  test.set  <- setdiff(rownames(K), train.set)
  
  # Setting the focal trait to NA for the test set:
  datalist$data.bm[datalist$data.bm$G %in% test.set, ncol(datalist$data.bm)] <- NA
  datalist$data.bm <- droplevels(datalist$data.bm)
  
  datalist$data.real[datalist$data.real$G %in% test.set, ncol(datalist$data.real)] <- NA
  datalist$data.real <- droplevels(datalist$data.real)
  
  # Subsetting pred.target:
  names(datalist$pred.target) <- rownames(K)
  datalist$pred.target <- datalist$pred.target[test.set]
  
  # Adding the training and test set vectors to datalist:
  datalist$test.set <- test.set
  datalist$train.set <- train.set
  
  # Scaling/centering training and test data separately:
  # Benchmark data:
  data.bm.test <- datalist$data.bm[which(is.na(datalist$data.bm$Y)),]
  data.bm.train <- datalist$data.bm[which(!is.na(datalist$data.bm$Y)),]
  
  data.bm.test[, 2:ncol(data.bm.test)] <- sapply(data.bm.test[, 2:ncol(data.bm.test)], scale)
  data.bm.train[, 2:ncol(data.bm.train)] <- sapply(data.bm.train[, 2:ncol(data.bm.train)], scale)
  
  datalist$data.bm <- rbind(data.bm.test, data.bm.train)
  
  # Real data:
  data.real.test <- datalist$data.real[which(is.na(datalist$data.real$Y)),]
  data.real.train <- datalist$data.real[which(!is.na(datalist$data.real$Y)),]
  
  data.real.test[, 2:ncol(data.real.test)] <- sapply(data.real.test[, 2:ncol(data.real.test)], scale)
  data.real.train[, 2:ncol(data.real.train)] <- sapply(data.real.train[, 2:ncol(data.real.train)], scale)
  
  datalist$data.real <- rbind(data.real.test, data.real.train)
  
  # Overwriting dataset:
  list.save(datalist, file = sprintf("timing/datasets/timing_p%d.RData", p))
  
  cat(sprintf("Generation of dataset for p = %s done!\n\n", p))
  toc()
}

end <- Sys.time()
end - start
