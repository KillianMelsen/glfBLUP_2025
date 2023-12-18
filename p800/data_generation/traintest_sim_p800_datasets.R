# For execution on WSL using oneMKL:
# Add `export MKL_NUM_THREADS=3` and `export MKL_DYNAMIC=FALSE` to ~/.profile (assuming a system with 20 threads)
# 5 parallel processes using foreach * 3 MKL threads = 15 out of 20 threads total.

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
    
    cl <- parallel::makeCluster(5, outfile = sprintf("logs/traintest_p800_h2s%s_comm%s.txt", h2s, comm))
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
        datalist <- list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData",
                                             h2y, comm, h2s, sim))
        
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
        list.save(datalist, file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData",
                                           h2y, comm, h2s, sim))
        
      }
    }
    toc()
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    combi <- combi + 1
  }
}


