# This script edits the datafiles generated using the script at
# `p800_lowrank/data_generation/generate_sim_p800_lowrank_datasets.R`. It selects
# the training and test sets and scales the data separately for training and test set.
#
# To check reproducibility, we can compare what training and test genotypes are
# selected for a specific replicate of any of the combinations of secondary
# feature heritability, communality, and focal trait heritability.
# Suppose we want to check replicate 96 for `h2s = "05"`, `comm = "05"`, and
# `h2y = "03"`. The first step is to set these values manually instead of relying
# on the different loops. Then we set `n.sim <- 96` instead of setting it to 100
# on line 69. Then, we comment out the code on lines 119/120 to prevent overwriting
# the datafile already saved to disk. We also comment out lines 77/78 that would load
# in the datafile already on disk. Finally we comment out all code from line 84 to 116.
# 
# We can now set the seed before running the inner loop. This will generate all the
# random training and test set divisions without actually modifying any of the datafiles
# already on disk. We can now compare the training and test sets that were sampled for
# replicate 96 to the training and test set according to the file on disk:
#
# all.equal(train.set,
#           rlist::list.load(sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))$train.set)
#
# [1] TRUE
#
# all.equal(test.set,
#           rlist::list.load(sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))$test.set)
#
# [1] TRUE

wd <- getwd()
setwd(wd)

# Libraries:
library(doParallel)
library(tictoc)

# Simulated genetic parameters:
h2.sec <- c("05", "07", "09")
comms <- c("02", "05", "08")
h2.foc <- c("01", "03", "05", "07", "09")

combis <- length(h2.sec) * length(comms)
combi <- 1

for (h2s in h2.sec) {
  for (comm in comms) {
    
    cl <- parallel::makeCluster(5, outfile = sprintf("logs/traintest_p800_lowrank_h2s%s_comm%s.txt", h2s, comm))
    doParallel::registerDoParallel(cl)
    
    cat(sprintf("Starting parallel train/test division of 5 x 100 simulated datasets (h2s = %s, comm = %s, combination %d/%d)...\n",
                h2s, comm, combi, combis))
    tic(sprintf("h2s = %s, comm = %s", h2s, comm))
    
    foreach::foreach(i = 1:length(h2.foc)) %dopar% {
      
      # Getting h2F chr value:
      h2y <- h2.foc[i]
      
      # Loading libraries:
      library(rlist)
      
      # Setting seed and working directory at start of 100 simulated datasets:
      set.seed(1997)
      setwd(wd)
      
      # Loading kinship and marker data, setting number of simulations and number of training genotypes:
      load("genotypes/K_sim.RData"); rm(M)
      n.sim <- 100
      n.train <- 300
      
      # Running train/test divisions on simulated datasets:
      sim <- 1
      for (sim in 1:n.sim) {
        
        # Loading raw simulated dataset:
        # datalist <- list.load(file = sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData",
        #                                      h2y, comm, h2s, sim))
        
        # Randomly choosing training and test genotypes:
        train.set <- sample(rownames(K), n.train)
        test.set  <- setdiff(rownames(K), train.set)
        
        # # Setting the focal trait to NA for the test set:
        # datalist$data.bm[datalist$data.bm$G %in% test.set, ncol(datalist$data.bm)] <- NA
        # datalist$data.bm <- droplevels(datalist$data.bm)
        # 
        # datalist$data.real[datalist$data.real$G %in% test.set, ncol(datalist$data.real)] <- NA
        # datalist$data.real <- droplevels(datalist$data.real)
        # 
        # # Subsetting pred.target:
        # names(datalist$pred.target) <- rownames(K)
        # datalist$pred.target <- datalist$pred.target[test.set]
        # 
        # # Adding the training and test set vectors to datalist:
        # datalist$test.set <- test.set
        # datalist$train.set <- train.set
        # 
        # # Scaling/centering training and test data separately:
        # # Benchmark data:
        # data.bm.test <- datalist$data.bm[which(is.na(datalist$data.bm$Y)),]
        # data.bm.train <- datalist$data.bm[which(!is.na(datalist$data.bm$Y)),]
        # 
        # data.bm.test[, 2:ncol(data.bm.test)] <- sapply(data.bm.test[, 2:ncol(data.bm.test)], scale)
        # data.bm.train[, 2:ncol(data.bm.train)] <- sapply(data.bm.train[, 2:ncol(data.bm.train)], scale)
        # 
        # datalist$data.bm <- rbind(data.bm.test, data.bm.train)
        # 
        # # Real data:
        # data.real.test <- datalist$data.real[which(is.na(datalist$data.real$Y)),]
        # data.real.train <- datalist$data.real[which(!is.na(datalist$data.real$Y)),]
        # 
        # data.real.test[, 2:ncol(data.real.test)] <- sapply(data.real.test[, 2:ncol(data.real.test)], scale)
        # data.real.train[, 2:ncol(data.real.train)] <- sapply(data.real.train[, 2:ncol(data.real.train)], scale)
        # 
        # datalist$data.real <- rbind(data.real.test, data.real.train)
        
        # Overwriting dataset:
        # list.save(datalist, file = sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData",
        #                                    h2y, comm, h2s, sim))
        
      }
    }
    toc()
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    combi <- combi + 1
  }
}


