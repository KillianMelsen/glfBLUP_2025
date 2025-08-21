# This script produces all of the simulated datasets with random residual
# structure (p800).
#
# The script is structured as follows:
# - There are two outer loops for secondary feature heritabilities ("05", "07", "09"),
#   and communality ("02", "05", "08").
# - There is then a parallel loop that divides the focal trait heritabilities
#   ("01", "03", "05", "07", "09") over 5 parallel workers.
# - The focal trait heritability for the datasets produced by each worker is assigned
#   on line 61.
# - The seed is then set within each worker and the final inner loop starting on
#   line 86 produces the 100 datasets for the replications for each combination
#   of secondary feature heritability, communality, and focal trait heritability.
#
# Suppose we want to reproduce dataset 13 for `h2s = "05"`, `comm = "08"`, and
# `h2y = "07"`. The first step is to set these values manually instead of relying
# on the two outer loops (h2s and comm) and the h2y assignment on line 61. Then we
# set `n.sim <- 13` instead of setting it to 100 on line 72. We then have to comment
# out the code on line 104/105 to avoid overwriting the dataset already saved to disk.
# We can then set the seed and run the inner loop. Finally, we can compare the
# `datalist` object we end up with to the object already saved to disk.
# Note that we can't directly compare the actual datasets because they still have
# to be scaled for training and test set using another script. Comparing the simulated
# genetic values will suffice to ensure reproducibility.
#
# all.equal(rlist::list.load(sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))$G,
#           datalist$G)
#
# [1] TRUE

start <- Sys.time()
wd <- getwd()

# Libraries:
library(doParallel)
library(tictoc)

# Simulated genetic parameters:
h2.sec <- c("05", "07", "09")
comms <- c("02", "05", "08")
h2.foc <- c("01", "03", "05", "07", "09")

combis <- length(h2.sec) * length(comms)
combi <- 1

h2s <- h2.sec[1]
comm <- comms[1]
for (h2s in h2.sec) {
  for (comm in comms) {
    
    cl <- parallel::makeCluster(5, outfile = sprintf("logs/generate_p800_h2s%s_comm%s.txt", h2s, comm))
    doParallel::registerDoParallel(cl)
    
    cat(sprintf("Starting parallel generation of 5 x 100 simulated datasets (h2s = %s, comm = %s, combination %d/%d)...\n",
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
      
      # Loading kinship and marker data, setting number of simulations and replicates:
      load("genotypes/K_sim.RData"); rm(M)
      n.sim <- 100
      r <- 2
      
      # Determining numeric values for the genetic parameters:
      comm.num <- as.numeric(comm) / 10
      
      sg2.s <- as.numeric(h2s) / 10
      se2.s <- 1 - sg2.s
      
      sg2.y <- as.numeric(h2y) / 10
      se2.y <- 1 - sg2.y
      
      # Running simulations:
      sim <- 1
      for (sim in 1:n.sim) {
        
        datalist <- glfBLUP::GFAsim(K = K, r = r, n.LSF = 4, n.LNF = 4,
                                    LSF.rho = 0, LNF.rho = 0, LSNF.rho = 0,
                                    S.per.LF = 100, Y.psi = 1 - comm.num,
                                    L.min = 0.3, L.max = 0.8,
                                    S.sg2 = sg2.s, S.se2 = se2.s,
                                    Y.sg2 = sg2.y, Y.se2 = se2.y,
                                    resCors = TRUE)
        
        datalist$simParams <- list(r = r, n.LSF = 4, n.LNF = 4,
                                   LSF.rho = 0, LNF.rho = 0, LSNF.rho = 0,
                                   S.per.LF = 100, Y.psi = 1 - comm.num,
                                   L.min = 0.3, L.max = 0.8,
                                   S.sg2 = sg2.s, S.se2 = se2.s,
                                   Y.sg2 = sg2.y, Y.se2 = se2.y,
                                   resCors = TRUE)
          
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

end <- Sys.time()
end - start


